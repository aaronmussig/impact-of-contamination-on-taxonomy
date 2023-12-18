import os

import dendropy
import luigi
from luigi import LocalTarget

from workflow.config import DIR_OUT_FASTTREE_FULL_TREE_AND_FAILED_NO_TRIM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_full_tree_and_failed_no_trim.b_run_fasttree import FastTreeFullTreeAndFailedNoTrimRunFastTree
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeFullTreeAndFailedNoTrimRoot(LuigiTask):
    replicate_id = luigi.IntParameter()

    @property
    def root_dir(self):
        if self.replicate_id == 0:
            return DIR_OUT_FASTTREE_FULL_TREE_AND_FAILED_NO_TRIM
        else:
            return os.path.join(DIR_OUT_FASTTREE_FULL_TREE_AND_FAILED_NO_TRIM, 'fasttree_replicates')

    def requires(self):
        return {
            'tree': FastTreeFullTreeAndFailedNoTrimRunFastTree(replicate_id=self.replicate_id),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'fasttree_undecorated_r{self.replicate_id}_rooted.tree'))

    def run(self):
        log(f'FastTreeFullTreeAndFailedNoTrimRoot(replicate={self.replicate_id})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading tree')
        tree = dendropy.Tree.get_from_path(self.input()['tree'].path, schema='newick', preserve_underscores=True)

        outgroup_gids = set(df_meta[df_meta['phylum'] == 'p__Spirochaetota'].index)
        outgroup_gids_with_prefix = {f'TEST_{x}' for x in outgroup_gids}

        outgroup_labels = {x.label for x in tree.taxon_namespace}.intersection(
            outgroup_gids.union(outgroup_gids_with_prefix))
        mrca_node = tree.mrca(taxon_labels=outgroup_labels)
        log(f'Found {len(mrca_node.leaf_nodes()):,} outgroup labels')

        tree.reroot_at_edge(mrca_node.edge, update_bipartitions=False)

        log(f'Writing to: {self.output().path}')
        tree.write(path=self.output().path, schema='newick', suppress_rooting=True, unquoted_underscores=True)
        return
