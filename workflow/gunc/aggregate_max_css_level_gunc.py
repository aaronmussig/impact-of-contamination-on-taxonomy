import multiprocessing as mp
import os

import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_GUNC, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc.run_gunc_on_gtdb import RunGuncOnGtdbR95UsingR207
from workflow.model.luigi import LocalTargetHdf5, LuigiTask
from workflow.util.gunc import MAX_CSS_LEVEL_DTYPE, get_max_css_path
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def worker(gid):
    root = get_gid_root(gid)
    path = get_max_css_path(os.path.join(root, 'gunc_r95'))
    df = pd.read_csv(path, sep='\t', dtype=MAX_CSS_LEVEL_DTYPE)
    if len(df) == 0:
        raise Exception(f'No data for {gid}')
    return df


class AggregateMaxCssLevelGtdbR95(LuigiTask):
    """Aggregates GUNC GTDB R95 MaxCSS into a single file."""

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'gunc': RunGuncOnGtdbR95UsingR207()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_GUNC, 'gunc_gtdb_95_max_css_level.h5'))

    def run(self):
        log('Aggregating max CSS level from GUNC directories (GTDB)', title=True)
        self.make_output_dirs()

        # Load the gids to process
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        gids = sorted(set(df_meta.index))
        log(f'{len(gids):,} genomes in GTDB R207')

        # Read each diamond file
        log('Reading GUNC all level files')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            results = list(tqdm(pool.imap_unordered(worker, gids), total=len(gids)))

        # Merge dataframe
        log('Concatenating data frames')
        df = pd.concat(results, ignore_index=True)
        df.rename(columns={"genome": "gid"}, inplace=True)

        self.save_hdf(df, index='gid')
