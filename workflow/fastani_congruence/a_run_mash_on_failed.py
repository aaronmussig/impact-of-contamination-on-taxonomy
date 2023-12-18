import multiprocessing as mp
import os
import tempfile

import luigi
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE, MASH_MIN_DIST, DEBUG
from workflow.external.mash_db import GtdbR207MashDb
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.genome_pct_congruence_contigs_removed import GenomePctCongruenceContigsRemoved
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.mash import SketchFile, DistanceFile
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class FastAniCongruenceRunMash(LuigiTask):
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'mash': GtdbR207MashDb(mash_k=self.mash_k, mash_s=self.mash_s),
            'contigs_removed': GenomePctCongruenceContigsRemoved(congruence=self.congruence, target_pct=self.target_pct),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(
            DIR_OUT_FASTANI_CONGRUENCE,
            f'mash_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}.h5'
        ))

    def run(self):
        log(f'FastAniCongruenceRunMash(p={self.target_pct}, c={self.congruence})', title=True)
        self.make_output_dirs()

        log('Loading contigs removed at congruence value')
        df_contigs_removed = self.input()['contigs_removed'].maybe_read_cached()

        log('Maybe caching mash sketch to disk')
        path_mash_sketch = self.input()['mash'].path_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].maybe_read_cached()

        log('Getting the pruned genome file for all failed')
        d_gid_to_pct_and_contigs_removed = dict()
        with tempfile.TemporaryDirectory() as tmp_dir:

            log('Creating queue')
            queue = list()
            for gid, row in tqdm(df_css.iterrows(), total=len(df_css)):
                if DEBUG and len(queue) > 10:
                    break

                cur_contigs_removed = frozenset(df_contigs_removed.loc[gid, 'contigs_removed'].split('|'))
                cur_pct_removed = float(df_contigs_removed.loc[gid, 'pct_removed'])

                queue.append((
                    gid,
                    cur_contigs_removed,
                    cur_pct_removed,
                    tmp_dir,
                ))

            log('Processing queue')
            with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
                for result in list(
                        tqdm(pool.imap_unordered(write_genome_fna_worker, queue), total=len(queue), smoothing=0.05)):
                    d_gid_to_pct_and_contigs_removed[result[0]] = (result[1], result[2])

            log('Running Mash')
            df = mash_worker(sorted(d_gid_to_pct_and_contigs_removed.keys()), path_mash_sketch, tmp_dir, self.mash_k,
                             self.mash_s)

            if not DEBUG:
                log('Writing to disk')
                self.save_hdf(df)

        return


def write_genome_fna_worker(job):
    gid, contigs, pct_removed, tmp_dir = job

    genome = Genome(gid)

    path_fna_out = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}.fna')
    os.makedirs(os.path.dirname(path_fna_out), exist_ok=True)

    # Write the modified FNA file to disk
    with open(path_fna_out, 'w') as f:
        for contig, seq in genome.d_fna.items():
            if contig not in contigs:
                f.write(f'>{contig}\n{seq}\n')

    return gid, pct_removed, contigs


def mash_worker(gids, path_r_sketch, tmp_dir, mash_k, mash_s):
    path_q_sketch = os.path.join(tmp_dir, 'query.msh')

    # Create the query sketch file
    d_gids = {x: os.path.join(get_gid_root(x, root=tmp_dir), f'{x}.fna') for x in gids}
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
