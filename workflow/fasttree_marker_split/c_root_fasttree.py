import os

import dendropy
import luigi
from luigi import LocalTarget

from workflow.config import DIR_OUT_FASTTREE_MARKER_SPLIT
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_marker_split.b_run_fasttree import FastTreeMarkerSplitRunFastTree
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeMarkerSplitRootFastTree(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_MARKER_SPLIT, f'pct_{self.target_pct}')

    def requires(self):
        return {
            'tree': FastTreeMarkerSplitRunFastTree(target_pct=self.target_pct),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'fasttree_undecorated_rooted.tree'))

    def run(self):
        log(f'FastTreeMarkerSplitRootFastTree(pct={self.target_pct})', title=True)
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
