import os

import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTTREE_MARKER_SPLIT, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_marker_split.e_decorate_fasttree import FastTreeMarkerSplitDecorateFastTree
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastTreeMarkerSplitJumpsBetweenHalves(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_MARKER_SPLIT, f'pct_{self.target_pct}')

    def requires(self):
        return {
            'split_tree': FastTreeMarkerSplitDecorateFastTree(target_pct=self.target_pct),
            'meta': GtdbMetadataR207(),
            'css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, 'num_jumps_between_halves.h5'))

    def run(self):
        log(f'FastTreeMarkerSplitJumpsBetweenHalves(p={self.target_pct})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_domain = df_meta['domain'].to_dict()

        log('Loading CSS')
        df_css = self.input()['css'].maybe_read_cached()
        fail_gids = {x for x in df_css.index if d_gid_to_domain[x] == 'd__Bacteria'}
        log(f'Found {len(fail_gids):,} failed genomes')

        if DEBUG:
            fail_gids = set(list(fail_gids)[:100])

        log('Loading split tree')
        split_tree = self.input()['split_tree'].read()

        log('Encoding split tree')
        split_tree.encode_bipartitions()

        log('Calculating jumps in split tree')
        d_gid_to_jumps = calculate_jumps_between_nodes(split_tree, fail_gids)

        log('Creating dataframe')
        rows = list()
        for gid, n_jumps in d_gid_to_jumps.items():
            rows.append({'gid': gid, 'n_jumps': n_jumps})
        df = pd.DataFrame(rows)

        if not DEBUG:
            log('Saving DF')
            self.save_hdf(df, 'gid')


def go_up_until_node_found(start_node, end_node):
    n = 0
    cur_node = start_node
    while cur_node is not None:
        if cur_node == end_node:
            return n
        cur_node = cur_node.parent_node
        n += 1
    raise Exception('Unable to find node')


def calculate_jumps_between_nodes(tree, gids):
    out = dict()

    gid_to_node = {x.taxon.label: x for x in tree.leaf_node_iter()}

    for gid in tqdm(gids):
        label_keep = f'TEST_{gid}'
        label_chim = f'TEST_{gid}_C'
        mrca = tree.mrca(taxon_labels=[label_keep, label_chim])

        n_jumps_keep = go_up_until_node_found(gid_to_node[label_keep], mrca)
        n_jumps_chim = go_up_until_node_found(gid_to_node[label_chim], mrca)

        out[gid] = n_jumps_keep + n_jumps_chim - 1

    return out
