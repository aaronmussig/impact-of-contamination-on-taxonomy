import multiprocessing as mp
import os
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from workflow.config import DIR_OUT_MARKER_REP, R207_MARKERS, DIR_OUT_MARKER_FAIL, \
    DEBUG, DIR_OUT_MARKER
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_fail_markers import GetFailMarkers
from workflow.marker.get_rep_markers import GetRepMarkers
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


def load_markers_from_file(path):
    out = defaultdict(lambda: 0)
    for record in SeqIO.parse(path, 'fasta'):
        gid, gene_id, num = record.id.split('|')
        out[gid] += 1
    return out


def process_marker_worker(job):
    marker, rep_path, fail_path = job
    d_results_rep = load_markers_from_file(rep_path)
    d_results_fail = load_markers_from_file(fail_path)
    return marker, {**d_results_rep, **d_results_fail}


class CalculateCntOfMarkersFailed(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            '_fail_markers': GetFailMarkers(),
            '_rep_markers': GetRepMarkers(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'count_of_marker_genes_for_failed_and_reps.h5'))

    def run(self):
        log('Calculating marker count for failed genomes', title=True)
        self.make_output_dirs()

        queue = list()
        for marker in R207_MARKERS:
            rep_path = os.path.join(DIR_OUT_MARKER_REP, f'{marker}.faa')
            fail_path = os.path.join(DIR_OUT_MARKER_FAIL, f'{marker}.faa')
            queue.append((marker, rep_path, fail_path))

        log('Processing queue')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(process_marker_worker, queue), total=len(queue)))

        log('Concatenating results')
        d_gid_to_marker_count = defaultdict(lambda: defaultdict(lambda: 0))
        for marker, d_results in tqdm(results):
            for gid, cnt in d_results.items():
                d_gid_to_marker_count[gid][marker] = cnt

        log('Creating dataframe rows')
        rows = list()
        for gid, d_marker_count in tqdm(d_gid_to_marker_count.items()):
            cur_row = {'gid': gid}
            for marker in R207_MARKERS:
                cur_row[marker] = d_marker_count[marker]
            rows.append(cur_row)

        log(f'Creating dataframe for {len(rows):,} rows')
        df = pd.DataFrame(rows)

        if not DEBUG:
            self.save_hdf(df, index='gid')
        return
