import multiprocessing as mp
import os
import tempfile

import luigi
from tqdm import tqdm

from workflow.config import MASH_MIN_DIST, DEBUG, DIR_OUT_FASTANI_CONTIG_SPLIT_RANDOM
from workflow.external.mash_db import GtdbR207MashDb
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.mash import SketchFile, DistanceFile
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class FastAniContigSplitRandomRunMash(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)
    replicate_id = luigi.IntParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTANI_CONTIG_SPLIT_RANDOM, f'r_{self.replicate_id}')

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'mash': GtdbR207MashDb(mash_k=self.mash_k, mash_s=self.mash_s),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(
            self.root_dir, f'mash_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}_r{self.replicate_id}.h5'
        ))

    def run(self):
        log(f'FastAniContigSplitRandomRunMash(p={self.target_pct}, r={self.replicate_id})', title=True)
        self.make_output_dirs()

        log('Maybe caching mash sketch to disk')
        path_mash_sketch = self.input()['mash'].path_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].maybe_read_cached()

        log('Getting the pruned genome file for all failed')
        with tempfile.TemporaryDirectory() as tmp_dir:

            log('Creating queue')
            queue = list()
            for i, (gid, row) in tqdm(enumerate(df_css.iterrows()), total=len(df_css)):
                if DEBUG and len(queue) > 10:
                    break

                queue.append((
                    gid,
                    self.target_pct,
                    i + self.replicate_id,
                    tmp_dir,
                ))

            log('Processing queue')
            with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
                results = list(
                    tqdm(pool.imap_unordered(write_genome_fna_worker, queue), total=len(queue), smoothing=0.05))

            log('Running Mash')
            df = mash_worker(results, path_mash_sketch, tmp_dir, self.mash_k, self.mash_s)

            if not DEBUG:
                log('Writing to disk')
                self.save_hdf(df)

        return


def write_genome_fna_worker(job):
    gid, pct_to_remove, seed, tmp_dir = job

    # Load the genome and split it into the target percent
    genome = Genome(gid)
    d_contig_to_seq_keep, d_contig_to_seq_disc = genome.split_fna_by_pct_random(pct=pct_to_remove, seed=seed)

    path_fna_keep = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}.fna')
    os.makedirs(os.path.dirname(path_fna_keep), exist_ok=True)

    path_fna_disc = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}_C.fna')
    os.makedirs(os.path.dirname(path_fna_disc), exist_ok=True)

    # Write the modified FNA files to disk
    with open(path_fna_keep, 'w') as f:
        for contig, seq in d_contig_to_seq_keep.items():
            f.write(f'>{contig}\n{seq}\n')

    with open(path_fna_disc, 'w') as f:
        for contig, seq in d_contig_to_seq_disc.items():
            f.write(f'>{contig}\n{seq}\n')

    return gid, path_fna_keep, path_fna_disc


def mash_worker(split_results, path_r_sketch, tmp_dir, mash_k, mash_s):
    path_q_sketch = os.path.join(tmp_dir, 'query.msh')

    # Create the query sketch file
    d_gids = dict()
    for result in split_results:
        d_gids[result[0]] = result[1]
        d_gids[f'{result[0]}_C'] = result[2]

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
    df_mash_data['query'] = df_mash_data['query'].apply(lambda x: os.path.basename(x).replace('.fna', ''))
    df_mash_data['ref'] = df_mash_data['ref'].apply(lambda x: os.path.basename(x).replace('.fna', ''))

    log('Sorting values')
    df_mash_data.sort_values(by=['query', 'dist', 'p_val'],
                             ascending=[True, True, True],
                             inplace=True,
                             ignore_index=True)

    return df_mash_data
