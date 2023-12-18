import os

import dendropy
import luigi
import pandas as pd
from phylodm import PhyloDM
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTTREE_MARKER_SPLIT_TRUE_CASE, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_marker_split_true_case.b_run_fasttree import FastTreeMarkerSplitTrueCaseRunFastTree
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.taxonomy import calculate_taxonomic_novelty


class FastTreeMarkerSplitTrueCaseCalculatePatristicDistance(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)
    replicate_id = luigi.IntParameter()

    @property
    def root_dir(self):
        return os.path.join(
            DIR_OUT_FASTTREE_MARKER_SPLIT_TRUE_CASE,
            f'pct_{self.target_pct}',
            f'replicate_{self.replicate_id}'
        )

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'tree': FastTreeMarkerSplitTrueCaseRunFastTree(target_pct=self.target_pct, replicate_id=self.replicate_id),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, f'patristic_distances.h5'))

    def run(self):
        log(f'FastTreeMarkerSplitTrueCaseCalculatePatristicDistance(target_pct={self.target_pct}, r={self.replicate_id})',
            title=True)
        self.make_output_dirs()

        log('Loading Tree')
        tree_path = self.input()['tree'].path

        log('Calculating Patristic Distance')
        pdm = generate_pdm(tree_path)
        all_gids = {x.replace('TEST_', '').replace('_C', '') for x in pdm.taxa()}
        fail_gids = {x.replace('TEST_', '').replace('_C', '') for x in pdm.taxa() if x.startswith('TEST_')}
        log(f'Found {len(fail_gids):,} test genome ids')

        log('Computing pairwise distances between each half')
        d_gid_to_pd = dict()
        tree_length = pdm.length()
        for gid in tqdm(sorted(fail_gids)):
            d_gid_to_pd[gid] = pdm.distance(f'TEST_{gid}', f'TEST_{gid}_C', norm=False) / tree_length

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_tax = df_meta[df_meta.index.isin(all_gids)]['gtdb_taxonomy'].to_dict()

        log('Computing taxonomic novelty for each genome')
        d_gid_to_tax_novelty, d_tax_novelty_to_gid = calculate_taxonomic_novelty(d_gid_to_tax)

        log('Recording values')
        rows = list()
        for gid in fail_gids:
            tax_novelty = d_gid_to_tax_novelty[gid]
            dist = d_gid_to_pd[gid]
            rows.append({
                'gid': gid,
                'tax_novelty': tax_novelty,
                'pd': dist,
            })
        df = pd.DataFrame(rows)

        if not DEBUG:
            log('Saving dataframe')
            self.save_hdf(df, index='gid')

        return


def generate_pdm(path_tree: str):
    tree = dendropy.Tree.get_from_path(path_tree, schema='newick', preserve_underscores=True)

    # if DEBUG:
    #     log('TESTING...')
    #     for taxon in tree.taxon_namespace:
    #         taxon.label = taxon.label[3:]
    #
    #     tree.taxon_namespace[0].label = f'TEST_{tree.taxon_namespace[5].label}_C'
    #     tree.taxon_namespace[5].label = f'TEST_{tree.taxon_namespace[5].label}'
    #
    #     tree.taxon_namespace[1].label = f'TEST_{tree.taxon_namespace[10].label}_C'
    #     tree.taxon_namespace[10].label = f'TEST_{tree.taxon_namespace[10].label}'

    n_updated = 0
    for edge in tree.postorder_edge_iter():
        if edge.length is None and not tree.seed_node.edge == edge:  # TODO: Make sure this is not the seed node edge
            edge.length = 0.0
            n_updated += 1
    log('Updated {} edges with missing lengths'.format(n_updated))

    log('Creating PDM')
    pdm = PhyloDM.load_from_dendropy(tree)

    log('Creating row vector')
    pdm.compute_row_vec()
    return pdm
