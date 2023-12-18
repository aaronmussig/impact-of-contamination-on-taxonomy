import multiprocessing as mp
import os
import tempfile

import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_MARKER
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.aggregate_max_css_level_merged_pass import AggregateMaxCssLevelMergedAll
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def worker(job):
    gid, max_css, source = job
    source = GuncRefDb.GTDB if source == 'gtdb' else GuncRefDb.PRO

    g = Genome(gid)
    d_marker_congruence = g.get_marker_congruence(max_css, source)

    out = list()
    for marker, congruence in d_marker_congruence.items():
        out.append({
            'gid': gid,
            'marker': marker,
            'congruence': congruence,
        })
    return out




class MarkerCongruenceFail(LuigiTask):
    """Contains the marker, genome, and congruence (pct correct) for all genomes"""

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'marker_congruence_all.h5'))

    def run(self):
        log('Getting congruence of all markers', title=True)
        self.make_output_dirs()

        df = self.input()['max_css'].maybe_read_cached()
        print(df.shape)

        queue = list()
        log('Creating queue')
        for gid, row in tqdm(df[['taxonomic_level', 'source']].iterrows(), total=len(df)):
            queue.append((gid, row['taxonomic_level'], row['source']))

            if DEBUG and len(queue) > 10:
                break

        log('Processing queue')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log('Converting to dataframe')
        rows = list()
        [rows.extend(x) for x in results]

        df = pd.DataFrame(rows)
        df.sort_values(by=['gid', 'marker'], inplace=True, ignore_index=True)

        log('Saving dataframe')
        if not DEBUG:
            self.save_hdf(df)

        return
