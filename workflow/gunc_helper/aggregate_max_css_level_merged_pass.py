import os

import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_GUNC, DEBUG
from workflow.gunc_helper.aggregate_max_css_level_gunc import AggregateMaxCssLevelGtdbR95
from workflow.gunc_helper.aggregate_max_css_level_progenomes import AggregateMaxCssLevelProGenomes
from workflow.model.luigi import LocalTargetHdf5, LuigiTask
from workflow.util.log import log


def merge_dfs(gtdb, pro):
    rows = list()

    common_gids = set(gtdb.index).intersection(set(pro.index))
    log(f'Found {len(common_gids):,} common gids')
    for gid in tqdm(common_gids):
        gtdb_row = gtdb.loc[gid]
        pro_row = pro.loc[gid]

        if gtdb_row.clade_separation_score >= pro_row.clade_separation_score:
            row = gtdb_row
        else:
            row = pro_row

        rows.append(row)

    return pd.DataFrame(rows)


class AggregateMaxCssLevelMergedAll(LuigiTask):
    """Unite both the GUNC and ProGenome MAX CSS into one file. Taking the MAX CSS."""

    def requires(self):
        return {
            'gunc_gtdb': AggregateMaxCssLevelGtdbR95(),
            'gunc_pro': AggregateMaxCssLevelProGenomes()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_GUNC, 'gunc_merged_max_css_level_all.h5'))

    def run(self):
        log('Uniting GUNC and ProGenome MAX CSS files (pass)', title=True)
        self.make_output_dirs()

        # Load the data frames
        df_gtdb = self.input()['gunc_gtdb'].read_cached() if DEBUG else self.input()['gunc_gtdb'].read()
        log(f'Found {len(df_gtdb):,} GTDB genomes')
        df_pro = self.input()['gunc_pro'].read_cached() if DEBUG else self.input()['gunc_pro'].read()
        log(f'Found {len(df_pro):,} ProGenomes genomes')

        # Include the soruce
        df_gtdb['source'] = 'gtdb'
        df_pro['source'] = 'progenomes'

        # Take the result with the highest CSS
        df_merged = merge_dfs(df_gtdb, df_pro)
        log(f'Found {len(df_merged):,} genomes (union)')

        n_gtdb = sum(df_merged['source'] == 'gtdb')
        n_pro = sum(df_merged['source'] == 'progenomes')
        log(f'{n_gtdb:,} were GTDB, {n_pro:,} were ProGenomes')

        # Write the data frame
        df_merged.reset_index(inplace=True)
        df_merged.rename(columns={"index": "gid"}, inplace=True)

        if not DEBUG:
            self.save_hdf(df_merged, index='gid')
