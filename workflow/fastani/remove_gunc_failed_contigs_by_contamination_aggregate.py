import multiprocessing as mp
import os

import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTANI
from workflow.fastani.remove_gunc_failed_contigs_by_contamination import RemoveGuncFailedContigsByContamination
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def worker(gid):
    gid_root = get_gid_root(gid)

    path_fastani = os.path.join(gid_root, 'fastani_gunc_failed_by_contamination.h5')

    if DEBUG and not os.path.isfile(path_fastani):
        return None

    df = pd.read_hdf(path_fastani)
    df['gid'] = gid
    df = df[['gid', 'reference', 'pct', 'ani', 'af']]
    return df


class RemoveGuncFailedContigsByContaminationAggregate(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            '_remove_contigs': RemoveGuncFailedContigsByContamination(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_FASTANI, 'remove_gunc_failed_contigs_by_contamination.h5'))

    def run(self):
        log('Removing contigs ordered by GUNC contamination (for failed only) (aggregate)', title=True)
        self.make_output_dirs()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Creating queue')
        gids = sorted(df_css.index)
        log(f'Queue size: {len(gids):,}')

        log('Aggregating files')
        if DEBUG:
            results = [worker(x) for x in gids[:100]]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(worker, gids), total=len(gids)))

        log('Concatenating dataframe')
        df = pd.concat(results)

        log('Saving dataframe')
        if not DEBUG:
            self.save_hdf(df)
