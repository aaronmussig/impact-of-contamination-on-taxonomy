import os

import pandas as pd

from workflow.config import DIR_OUT_GUNC, DEBUG
from workflow.gunc.aggregate_max_css_level_gunc import AggregateMaxCssLevelGtdbR95
from workflow.gunc.aggregate_max_css_level_progenomes import AggregateMaxCssLevelProGenomes
from workflow.model.luigi import LocalTargetHdf5, LuigiTask
from workflow.util.log import log


def merge_dfs(gtdb, pro):
    rows = list()

    common_gids = set(gtdb.index).intersection(set(pro.index))
    log(f'Found {len(common_gids):,} common gids that failed')
    for gid in common_gids:
        gtdb_row = gtdb.loc[gid]
        pro_row = pro.loc[gid]

        if gtdb_row.clade_separation_score >= pro_row.clade_separation_score:
            row = gtdb_row
        else:
            row = pro_row

        rows.append(row)

    for gid, row in gtdb.loc[~gtdb.index.isin(common_gids)].iterrows():
        rows.append(row)

    for gid, row in pro.loc[~pro.index.isin(common_gids)].iterrows():
        rows.append(row)

    return pd.DataFrame(rows)


class AggregateMaxCssLevelMerged(LuigiTask):
    """Unite both the GUNC and ProGenome MAX CSS into one file. Taking the MAX CSS."""

    def requires(self):
        return {
            'gunc_gtdb': AggregateMaxCssLevelGtdbR95(),
            'gunc_pro': AggregateMaxCssLevelProGenomes()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_GUNC, 'gunc_merged_max_css_level.h5'))

    def run(self):
        log('Uniting GUNC and ProGenome MAX CSS files', title=True)
        self.make_output_dirs()

        # Load the data frames
        df_gtdb = self.input()['gunc_gtdb'].read_cached() if DEBUG else self.input()['gunc_gtdb'].read()
        log(f'Found {len(df_gtdb):,} GTDB genomes')
        df_pro = self.input()['gunc_pro'].read_cached() if DEBUG else self.input()['gunc_pro'].read()
        log(f'Found {len(df_pro):,} ProGenomes genomes')

        # Reduce to the subset that failed
        df_gtdb = df_gtdb[df_gtdb['pass.GUNC'] == False]
        log(f'Found {len(df_gtdb):,} GTDB genomes that failed')
        df_pro = df_pro[df_pro['pass.GUNC'] == False]
        log(f'Found {len(df_pro):,} ProGenomes genomes that failed')

        # Include the soruce
        df_gtdb['source'] = 'gtdb'
        df_pro['source'] = 'progenomes'

        # Take the result with the highest CSS
        df_merged = merge_dfs(df_gtdb, df_pro)
        log(f'Found {len(df_merged):,} genomes that failed (union)')

        n_gtdb = sum(df_merged['source'] == 'gtdb')
        n_pro = sum(df_merged['source'] == 'progenomes')
        log(f'{n_gtdb:,} were GTDB, {n_pro:,} were ProGenomes')

        # Write the data frame
        df_merged.reset_index(inplace=True)
        df_merged.rename(columns={"index": "gid"}, inplace=True)
        self.save_hdf(df_merged, index='gid')
