import os

import dendropy
import luigi
from luigi import LocalTarget

from workflow.config import DIR_OUT_FASTTREE_MARKER_SPLIT_RANDOM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_marker_split.random_b_run_fasttree import FastTreeMarkerSplitRandomRunFastTree
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeMarkerSplitRandomRootFastTree(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)
    replicate_id = luigi.IntParameter()

    @property
    def root_dir(self):
        return os.path.join(
            DIR_OUT_FASTTREE_MARKER_SPLIT_RANDOM,
            f'pct_{self.target_pct}',
            f'replicate_{self.replicate_id}'
        )

    def requires(self):
        return {
            'tree': FastTreeMarkerSplitRandomRunFastTree(target_pct=self.target_pct, replicate_id=self.replicate_id),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'fasttree_undecorated_rooted.tree'))

    def run(self):
        log(f'FastTreeMarkerSplitRandomRootFastTree(pct={self.target_pct}, r={self.replicate_id})', title=True)
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
