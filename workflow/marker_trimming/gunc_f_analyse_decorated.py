import os

import luigi
from luigi import LocalTarget

from workflow.config import DIR_OUT_MARKER_TRIMMING_GUNC
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_full_tree_and_failed_no_trim.e_decorate_fasttree import FastTreeFullTreeAndFailedNoTrimDecorate
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker_trimming.gunc_e_decorate_fasttree import MarkerTrimmingGuncDecorateFastTree
from workflow.method.check_phylorank_decorate_output import check_phylorank_decorate_output
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class MarkerTrimmingGuncDecorateAnalyseFastTree(LuigiTask):
    target_pct = luigi.FloatParameter()

    @property
    def root_dir(self):
        return DIR_OUT_MARKER_TRIMMING_GUNC

    def requires(self):
        return {
            'decorated': MarkerTrimmingGuncDecorateFastTree(target_pct=self.target_pct),
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'ref_tree': FastTreeFullTreeAndFailedNoTrimDecorate(),
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'p{self.target_pct}', 'analysis_results_tmp.tsv'))

    def run(self):
        log(f'MarkerTrimmingGuncDecorateAnalyseFastTree(pct={self.target_pct})', title=True)
        self.make_output_dirs()

        log('Loading BAC120 R207 reference tree')
        ref_tree = self.input()['ref_tree'].read()
        for leaf_node in ref_tree.leaf_node_iter():
            leaf_node.taxon.label = leaf_node.taxon.label[3:]

        log('Loading decorated tree')
        decorated_tree = self.input()['decorated'].read()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading decorated tree taxonomy')
        phylo_dir = os.path.dirname(self.input()['decorated'].path)
        path_decorated_tax = os.path.join(phylo_dir, 'decorated.tree-taxonomy')

        log('Comparing trees')
        results = check_phylorank_decorate_output(
            df_meta,
            ref_tree,
            decorated_tree,
            path_decorated_tax
        )

        return
