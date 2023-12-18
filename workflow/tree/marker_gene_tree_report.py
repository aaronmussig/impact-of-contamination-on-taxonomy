import os
import multiprocessing as mp
import os
import subprocess
import tempfile
from collections import defaultdict

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
from phylodm import PhyloDM
import dendropy


def process_marker(marker, df_meta):
    path_out = os.path.join(DIR_OUT_TREE_MARKER_GENE_TREE, marker)
    path_tree = os.path.join(path_out, f'{marker}.tree')

    tree = dendropy.Tree.get_from_path(path_tree, schema='newick', preserve_underscores=True)
    pdm = PhyloDM.load_from_dendropy(tree)

    dm = pdm.dm(norm=False)
    labels = pdm.taxa()

    dm_sorted = np.argsort(dm, axis=1)

    label_tax = [df_meta.loc[x[0:15], 'gtdb_taxonomy'] for x in labels]

    rep_gids = frozenset(df_meta[df_meta['gtdb_representative'] == 't'].index)

    for i, (label, tax) in enumerate(zip(labels, label_tax)):

        d_taxon_to_cnt = defaultdict(lambda: defaultdict(lambda: 0))
        n_gids_seen = 0

        for cur_idx in dm_sorted[i, :]:
            cur_gid = labels[cur_idx][0:15]
            cur_gid_is_rep = cur_gid in rep_gids

            if not cur_gid_is_rep:
                continue
            cur_taxon = label_tax[cur_idx]

            for rank in cur_taxon.split(';'):
                d_taxon_to_cnt[rank[0]][rank] += 1

            n_gids_seen += 1
            if n_gids_seen > 50:
                break



        print()
    return


class MarkerGeneTreeReport(LuigiTask):

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
        log('Creating gene trees from markers (report)', title=True)
        self.make_output_dirs()

        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        queue = list()
        for marker in tqdm(R207_MARKERS):

            if DEBUG and marker != 'PF01200.19':
                continue

            process_marker(marker, df_meta)

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
