import os
import multiprocessing as mp
import os
import subprocess
import tempfile

import h5py
import numpy as np
from Bio import SeqIO
from luigi import LocalTarget
from magna.util.disk import move_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DIR_OUT_SENTINEL, DIR_OUT_MARKER_REP, R207_MARKERS, DIR_OUT_MARKER_FAIL, \
    DIR_OUT_MARKER_AAI_FAIL_TO_REP, DEBUG, DIR_OUT_TREE_MARKER_GENE_TREE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_fail_markers import GetFailMarkers
from workflow.marker.get_rep_markers import GetRepMarkers
from workflow.model.luigi import LuigiTask
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.rq import rq_and_wait


def merge_marker_files(marker):
    path_fail = os.path.join(DIR_OUT_MARKER_FAIL, f'{marker}.faa')
    path_rep = os.path.join(DIR_OUT_MARKER_REP, f'{marker}.faa')

    d_fail = read_fasta(path_fail)
    d_rep = read_fasta(path_rep)
    return {**d_fail, **d_rep}

def mask_seqs(d_faa):

    arr = np.array([list(x) for x in d_faa.values()])

    mask = np.zeros((1, arr.shape[1]), dtype=bool)
    for i in range(arr.shape[1]):
        col_unq = np.unique(arr[:, 0])
        col_unq = {str(x) for x in col_unq}
        col_unq -= {'-'}

        if len(col_unq) == 1:
            mask[0, i] = True
            print('masked')

    return d_faa


def run_iqtree_on_marker(marker):

    path_out = os.path.join(DIR_OUT_TREE_MARKER_GENE_TREE, marker)
    os.makedirs(path_out, exist_ok=True)

    d_faa =merge_marker_files(marker)
    # d_faa_masked = mask_seqs(d_faa)

    path_faa_iqtree = os.path.join(path_out, f'{marker}.faa')
    with open(path_faa_iqtree, 'w') as f:
        for k, v in sorted(d_faa.items()):
            f.write(f'>{k}\n{v}\n')

    cmd = [
        'FastTreeMP',
        '-wag', path_faa_iqtree,
        '>', os.path.join(path_out, f'{marker}.tree')
    ]

    # print(' '.join(cmd))

    return


def create_tree(marker):
    path_out = os.path.join(DIR_OUT_TREE_MARKER_GENE_TREE, marker)
    os.makedirs(path_out, exist_ok=True)

    path_faa_iqtree = os.path.join(path_out, f'{marker}.faa')
    path_tree = os.path.join(path_out, f'{marker}.tree')
    path_log = os.path.join(path_out, f'{marker}.log')

    cmd = [
        'FastTreeMP',
        '-wag', path_faa_iqtree,
    ]

    with open(path_tree, 'w') as f_tree, open(path_log, 'w') as f_log:
        proc = subprocess.Popen(cmd, stdout=f_tree, stderr=f_log, encoding='utf-8')
        stdout, stderr = proc.communicate()

    return


class MarkerGeneTree(LuigiTask):

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
        log('Creating gene trees from markers', title=True)
        self.make_output_dirs()

        queue = list()
        for marker in tqdm(R207_MARKERS):
            queue.append((marker, ))
            # run_iqtree_on_marker(marker)

        # with mp.Pool(processes=min(mp.cpu_count(), 20)) as pool:
        #     list(tqdm(pool.imap_unordered(run_iqtree_on_marker, queue), total=len(queue)))
        #
        # if not DEBUG:
        #     self.write_sentinel()

        if not DEBUG:
            rq_and_wait(job_id=self.fqn, fn=create_tree, q_args=queue, queue_name=self.fqn)
            self.write_sentinel()

        return
