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
    DIR_OUT_MARKER_AAI_FAIL_TO_REP, DEBUG, R207_MARKER_LENGTHS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_fail_markers import GetFailMarkers
from workflow.marker.get_rep_markers import GetRepMarkers
from workflow.model.aai_fail_to_rep_h5 import read_aai_fail_to_rep_h5
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


def calculate_on_h5_path(marker, df_meta):
    path = os.path.join(DIR_OUT_MARKER_AAI_FAIL_TO_REP, f'{marker}.h5')
    marker_len = R207_MARKER_LENGTHS[marker]

    fail_keys, rep_keys, n_correct, n_aminos = read_aai_fail_to_rep_h5(path)
    rep_taxa = [df_meta.loc[x[0:15], 'gtdb_taxonomy'] for x in rep_keys]

    pct_correct = (n_correct / n_aminos) * 100
    pct_correct_sorted = np.argsort(pct_correct, axis=1)

    for i, fail_key in enumerate(fail_keys):
        fail_tax = df_meta.loc[fail_key[0:15], 'gtdb_taxonomy']
        out = list()

        for j in pct_correct_sorted[i, :]:
            rep_tax = rep_taxa[j]
            rep_key = rep_keys[j]

            if rep_key[0:15] == fail_key[0:15]:
                continue

            pct = pct_correct[i, j]

            if pct < 80:
                continue

            highest_agree_tax = ';'.join([f for f, r in zip(fail_tax.split(';'), rep_tax.split(';')) if f == r]).split(';')[-1][0]
            out.append((pct, highest_agree_tax, rep_tax))

        print()
    return


class CalculateAaiOfMarkersForFailedReport(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Calculating marker AAI for failed genomes', title=True)
        self.make_output_dirs()
        os.makedirs(DIR_OUT_MARKER_AAI_FAIL_TO_REP, exist_ok=True)

        log('Loading metadata')
        df_meta = self.input()['meta'].read_cached() if not DEBUG else self.input()['meta'].read()

        for marker in R207_MARKERS:
            if marker != 'PF04919.13':
                continue

            calculate_on_h5_path(marker, df_meta)

        return
