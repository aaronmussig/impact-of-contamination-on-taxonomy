import os
import random

import pandas as pd

from workflow.config import DEBUG, DIR_OUT_BOOTSTRAP
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


N_REPS = 10


class SelectGenomesForBootstrapping(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_BOOTSTRAP, 'select_genomes_for_bootstrapping.h5'))

    def run(self):
        log('Selecting genomes for bootstrapping', title=True)
        self.make_output_dirs()

        random.seed(42)

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        gids = sorted(df_meta.index)

        log('Loading MaxCSS merged file')
        df_max_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log(f'Selecting genomes from {len(gids):,} genomes @ {len(df_max_css):,} per sample ({N_REPS} repetitions)')
        rows = list()
        for n in range(N_REPS):
            rows.append({
                'rep': n,
                'genomes': '|'.join(random.choices(gids, k=len(df_max_css)))
            })

        log('Saving to disk')
        df = pd.DataFrame(rows)
        if not DEBUG:
            self.save_hdf(df, index='rep')
