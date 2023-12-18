import os
import multiprocessing as mp
import os
import tempfile

import h5py
import numpy as np
from Bio import SeqIO
from luigi import LocalTarget
from magna.util.disk import move_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DIR_OUT_SENTINEL, DIR_OUT_MARKER_REP, R207_MARKERS, DIR_OUT_MARKER_FAIL, \
    DIR_OUT_MARKER_AAI_FAIL_TO_REP, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_fail_markers import GetFailMarkers
from workflow.marker.get_rep_markers import GetRepMarkers
from workflow.model.luigi import LuigiTask
from workflow.util.log import log

AA = {'-': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8,
      'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16,
      'T': 17, 'V': 18, 'W': 19, 'X': 20, 'Y': 21}

MAX_HEAP_SIZE = 1000


def load_markers_from_file(path):
    # log(f'Loading markers from: {path}')
    keys = list()
    rows = list()
    for record in sorted(SeqIO.parse(path, 'fasta'), key=lambda x: x.id):
        keys.append(record.id)
        rows.append(str(record.seq))

    arr = np.zeros((len(rows), len(rows[0])), dtype=np.uint8)
    for i, row in enumerate(rows):
        for j, char in enumerate(row):
            arr[i, j] = AA[char]
    return tuple(keys), arr


def process_marker_worker(job):
    marker, rep_path, fail_path = job

    rep_keys, rep_arr = load_markers_from_file(rep_path)
    fail_keys, fail_arr = load_markers_from_file(fail_path)

    rep_arr_not_gap = rep_arr != 0
    fail_arr_not_gap = fail_arr != 0

    # log(f'Calculating AAI for: {marker}')

    marker_len = rep_arr.shape[1]
    if marker_len <= 255:
        dtype = np.uint8
    elif marker_len <= 65535:
        dtype = np.uint16
    else:
        raise Exception('????')

    # log(f'Creating output array: {dtype}')
    out_n_aminos = np.zeros((len(fail_keys), len(rep_keys)), dtype=dtype)
    out_n_correct = np.zeros((len(fail_keys), len(rep_keys)), dtype=dtype)

    for i, fail_key in enumerate(fail_keys):
        eq_values = fail_arr[i, :] == rep_arr
        cols_without_gap = fail_arr_not_gap[i, :] & rep_arr_not_gap

        n_aminos_correct = (eq_values & cols_without_gap).sum(axis=1)
        n_aminos_per_row = cols_without_gap.sum(axis=1)

        for j, rep_key in enumerate(rep_keys):
            out_n_aminos[i, j] = int(n_aminos_per_row[j])
            out_n_correct[i, j] = int(n_aminos_correct[j])

    path_out = os.path.join(DIR_OUT_MARKER_AAI_FAIL_TO_REP, f'{marker}.h5')
    with tempfile.TemporaryDirectory() as tmp_dir:
        path_df_tmp = os.path.join(tmp_dir, 'df.h5')

        # log(f'Saving data to: {path_df_tmp}')
        with h5py.File(path_df_tmp, 'w') as hf:
            hf.create_dataset('fail_keys', data=fail_keys, compression="gzip", compression_opts=9)
            hf.create_dataset('rep_keys', data=rep_keys, compression="gzip", compression_opts=9)
            hf.create_dataset('n_aminos', data=out_n_aminos, compression="gzip", compression_opts=9)
            hf.create_dataset('n_correct', data=out_n_correct, compression="gzip", compression_opts=9)

        if not DEBUG:
            log(f'Copying {get_file_size_fmt(path_df_tmp)} to: {path_out}')
            move_file(path_df_tmp, path_out, checksum=True)
            # log('Done.')

    return


class CalculateAaiOfMarkersForFailed(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            '_fail_markers': GetFailMarkers(),
            '_rep_markers': GetRepMarkers(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Calculating marker AAI for failed genomes', title=True)
        self.make_output_dirs()
        os.makedirs(DIR_OUT_MARKER_AAI_FAIL_TO_REP, exist_ok=True)

        queue = list()
        for marker in R207_MARKERS:
            rep_path = os.path.join(DIR_OUT_MARKER_REP, f'{marker}.faa')
            fail_path = os.path.join(DIR_OUT_MARKER_FAIL, f'{marker}.faa')
            queue.append((marker, rep_path, fail_path))

        with mp.Pool(processes=min(mp.cpu_count(), 20)) as pool:
            list(tqdm(pool.imap_unordered(process_marker_worker, queue), total=len(queue)))

        if not DEBUG:
            self.write_sentinel()

        return
