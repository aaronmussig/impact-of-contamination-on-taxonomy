import multiprocessing as mp
import os
import tempfile

import luigi
from tqdm import tqdm

from workflow.config import MASH_MIN_DIST, DEBUG, DIR_OUT_FASTANI_RANDOM
from workflow.external.mash_db import GtdbR207MashDb
from workflow.fastani_random.a_randomly_select_gids import FastAniRandomSelectGids
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.mash import SketchFile, DistanceFile
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class FastAniRandomRunMashOnRandomSelection(LuigiTask):
    target_pct = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    batch_id = luigi.IntParameter()

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'mash': GtdbR207MashDb(mash_k=self.mash_k, mash_s=self.mash_s),
            'gids': FastAniRandomSelectGids(target_pct=self.target_pct, batch_id=self.batch_id),
        }

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTANI_RANDOM, f'b{self.batch_id}__p{self.target_pct}')

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, 'b_mash_to_r207_reps.h5'))

    def run(self):
        log(f'FastAniRandomRunMashOnRandomSelection(b={self.batch_id}, p={self.target_pct})', title=True)
        self.make_output_dirs()

        log('(Maybe) caching the mash sketch file')
        path_mash_sketch = self.input()['mash'].path_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].maybe_read_cached()

        log('Loading randomly selected genomes')
        df_gids = self.input()['gids'].maybe_read_cached()

        log('Getting the pruned genome file for all failed')
        with tempfile.TemporaryDirectory() as tmp_dir:

            log('Creating queue')
            queue = list()
            for row in tqdm(df_gids.itertuples(), total=len(df_css)):
                if DEBUG and len(queue) > 10:
                    break

                uid = row.Index
                gid = row.gid
                contigs_removed = frozenset(row.contigs.split('|'))
                pct_removed = float(row.pct_removed)

                if pct_removed == 0.0 or len(contigs_removed) == 0:
                    continue

                queue.append((
                    uid,
                    gid,
                    contigs_removed,
                    pct_removed,
                    tmp_dir,
                ))

            log('Processing queue')
            mash_queue = list()
            with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
                for cur_uid, cur_gid, cur_path_fna_out in list(
                        tqdm(pool.imap_unordered(write_genome_fna_worker, queue), total=len(queue), smoothing=0.05)):
                    mash_queue.append((cur_uid, cur_gid, cur_path_fna_out))

            log('Running Mash')
            df = mash_worker(mash_queue, path_mash_sketch, tmp_dir, self.mash_k, self.mash_s)

            if not DEBUG:
                log('Writing to disk')
                self.save_hdf(df)

        return


def write_genome_fna_worker(job):
    uid, gid, contigs, pct_removed, tmp_dir = job

    genome = Genome(gid)

    path_fna_out = os.path.join(get_gid_root(gid, root=tmp_dir), f'{uid}_{gid}.fna')
    os.makedirs(os.path.dirname(path_fna_out), exist_ok=True)

    # Write the modified FNA file to disk
    with open(path_fna_out, 'w') as f:
        for contig, seq in genome.d_fna.items():
            if contig not in contigs:
                f.write(f'>{contig}\n{seq}\n')

    return uid, gid, path_fna_out


def mash_worker(mash_queue, path_r_sketch, tmp_dir, mash_k, mash_s):
    path_q_sketch = os.path.join(tmp_dir, 'query.msh')

    d_gids = dict()
    for cur_uid, cur_gid, cur_path_fna_out in mash_queue:
        d_gids[cur_uid] = cur_path_fna_out

    # Create the query sketch file
    SketchFile(genomes=d_gids, path=path_q_sketch, cpus=mp.cpu_count(), k=mash_k, s=mash_s)

    # Find the closest 100 representative genomes (by mash)
    mash = DistanceFile(
        qry_sketch=path_q_sketch,
        ref_sketch=path_r_sketch,
        cpus=mp.cpu_count(),
        max_d=MASH_MIN_DIST,
        mash_v=1.0
    )
    log('Getting mash data')
    df_mash_data = mash.get_data()
    df_mash_data['uid'] = df_mash_data['query'].apply(lambda x: get_uid_and_gid_from_tmp_path(x)[0])
    df_mash_data['query'] = df_mash_data['query'].apply(lambda x: get_uid_and_gid_from_tmp_path(x)[1])
    df_mash_data['ref'] = df_mash_data['ref'].apply(lambda x: os.path.basename(x).replace('.fna', ''))

    log('Sorting values')
    df_mash_data.sort_values(by=['uid', 'dist', 'p_val'],
                             ascending=[True, True, True],
                             inplace=True,
                             ignore_index=True)

    return df_mash_data


def get_uid_and_gid_from_tmp_path(path: str):
    name = os.path.basename(path)
    idx = name.index('_')
    uid = int(name[:idx])
    gid = name[idx + 1:-4]
    return uid, gid
