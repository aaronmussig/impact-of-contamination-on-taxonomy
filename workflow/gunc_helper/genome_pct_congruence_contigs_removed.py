import multiprocessing as mp
import os

import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_GUNC
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class GenomePctCongruenceContigsRemoved(LuigiTask):
    """GUNC-informed removal of contigs."""
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                DIR_OUT_GUNC,
                f'genome_p{self.target_pct}_c{self.congruence}_contigs_removed.h5'
            ))

    def run(self):
        log(f'GenomePctCongruenceContigsRemoved(p={self.target_pct}, c={self.congruence})', title=True)
        self.make_output_dirs()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].maybe_read_cached()

        log('Getting the pruned genome file for all failed')
        d_gid_to_pct_and_contigs_removed = dict()

        log('Creating queue')
        queue = list()
        for gid, row in tqdm(df_css.iterrows(), total=len(df_css)):
            queue.append((
                gid,
                row['source'],
                row['taxonomic_level'],
                self.congruence,
                self.target_pct,
            ))
            if DEBUG and len(queue) >= 10:
                break

        log('Processing queue')
        rows = list()
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            for result in list(tqdm(pool.imap_unordered(genome_fna_worker, queue), total=len(queue), smoothing=0.05)):
                d_gid_to_pct_and_contigs_removed[result[0]] = (result[1], result[2])
                rows.append({
                    'gid': result[0],
                    'pct_removed': result[1],
                    'contigs_removed': '|'.join(sorted(result[2])),
                })

        log('Creating output dataframe')
        df_out = pd.DataFrame(rows)

        log('Saving')
        if not DEBUG:
            self.save_hdf(df_out, index='gid')

        return


def genome_fna_worker(job):
    gid, source, max_css, congruence, target_pct = job

    if source == 'gtdb':
        source = GuncRefDb.GTDB
    elif source == 'progenomes':
        source = GuncRefDb.PRO
    else:
        raise ValueError(f'Unknown source: {source}')

    genome = Genome(gid)
    contigs, pct_removed = genome.get_gunc_contigs_where_removed_equals_x_pct_genome_removed(
        pct=target_pct,
        max_congruence=congruence,
        max_css=max_css,
        source=source
    )
    contigs = frozenset(contigs)

    return gid, pct_removed, contigs
