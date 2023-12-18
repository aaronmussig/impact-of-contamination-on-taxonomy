import os
import subprocess
import tempfile

import dendropy
import luigi
from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, DIR_OUT_FASTTREE_MARKER_ANALYSIS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_full_tree_non_reps.gunc_a_create_batches import FastTreeFullSetNonRepsCreateBatchesGunc
from workflow.fasttree_full_tree_non_reps.gunc_b_run_fasttree import FastTreeFullSetNonRepsCreateBatchesGuncFastTree
from workflow.fasttree_marker_analysis.gunc_b_run_fasttree import FastTreeMarkerAnalysisRunGuncFastTree
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeMarkerAnalysisRootGuncFastTree(LuigiTask):
    congruence = luigi.FloatParameter()
    target_pct = luigi.FloatParameter()

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS, f'c{str(self.congruence)}_p{str(self.target_pct)}')

    def requires(self):
        return {
            'tree': FastTreeMarkerAnalysisRunGuncFastTree(congruence=self.congruence, target_pct=self.target_pct),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, 'fasttree_undecorated_rooted.tree'))

    def run(self):
        log(f'FastTreeMarkerAnalysisRootGuncFastTree(congruence={self.congruence}) (pct={self.target_pct})',
            title=True)
        self.make_output_dirs()

        raise Exception("TODO!")

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading tree')
        tree = dendropy.Tree.get_from_path(self.input()['tree'].path, schema='newick', preserve_underscores=True)

        outgroup_gids = set(df_meta[df_meta['phylum'] == 'p__Spirochaetota'].index)
        outgroup_gids_with_prefix = {f'TEST_{x}' for x in outgroup_gids}

        outgroup_labels = {x.label for x in tree.taxon_namespace}.intersection(outgroup_gids.union(outgroup_gids_with_prefix))
        mrca_node = tree.mrca(taxon_labels=outgroup_labels)
        log(f'Found {len(mrca_node.leaf_nodes()):,} outgroup labels')

        tree.reroot_at_edge(mrca_node.edge, update_bipartitions=False)

        log(f'Writing to: {self.output().path}')
        tree.write(path=self.output().path, schema='newick', suppress_rooting=True, unquoted_underscores=True)
        return

