import multiprocessing as mp
import os
import random

import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTANI_RANDOM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastAniRandomSelectGids(LuigiTask):
    target_pct = luigi.FloatParameter()
    batch_id = luigi.IntParameter()

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTANI_RANDOM, f'b{self.batch_id}__p{self.target_pct}')

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, 'a_selected_genomes_and_contigs_removed.h5'))

    def run(self):
        log(f'FastAniRandomSelectGids(b={self.batch_id}, p={self.target_pct})', title=True)
        self.make_output_dirs()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].maybe_read_cached()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log(f'Getting a random selection of n={len(df_css):,} genomes')
        gids = random.choices(df_meta.index, k=len(df_css))
        queue = [(i, x, self.target_pct) for i, x in enumerate(sorted(gids))]
        queue = queue[:10] if DEBUG else queue

        log('Determining the contigs to remove from each genome')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            rows = list(tqdm(pool.imap_unordered(randomly_select_contigs_for_removal, queue), total=len(queue)))

        log('Creating dataframe')
        df = pd.DataFrame(rows)

        if not DEBUG:
            self.save_hdf(df, index='uid')

        return


def randomly_select_contigs_for_removal(job):
    job_id, gid, target_pct = job

    genome = Genome(gid)
    contigs, pct_removed = genome.get_random_contigs_where_removed_equals_x_pct_genome_removed(target_pct)

    return {
        'uid': job_id,
        'gid': gid,
        'pct_removed': pct_removed,
        'contigs': '|'.join(sorted(contigs)),
    }
