from collections import defaultdict
from typing import Tuple

import dendropy
import pandas as pd
from tqdm import tqdm

from workflow.fasttree_marker_analysis.gunc_c_decorate_fasttree import FastTreeMarkerAnalysisDecorateGuncFastTree
from workflow.util.tree import parse_label


def check_phylorank_decorate_output(
        df_meta: pd.DataFrame,
        ref_tree: dendropy.Tree,
        dec_tree: dendropy.Tree,
        path_decorated_tax: str,

):
    """
    This function takes a decorated input tree and checks the taxonomic
    assignment of prefixed genomes.
    """

    # Get the true taxonomy of each accession
    d_gid_to_true_taxonomy = df_meta['gtdb_taxonomy'].to_dict()

    # Get the decorated taxonomy of each accession
    d_gid_to_decorated_tax = _read_decorated_tax_file(path_decorated_tax)

    # Get the test gids
    test_gids = frozenset({x[5:] for x in d_gid_to_decorated_tax if x.startswith('TEST_')})

    # Iterate over all decorated taxonomy strings to check the TEST_ genomes
    wrong, correct = set(), set()
    for raw_gid, decorated_tax in tqdm(sorted(d_gid_to_decorated_tax.items())):
        if raw_gid.startswith('TEST_'):
            gid = raw_gid[5:]
        else:
            continue

        # Load the expected taxonomy string for this genome
        true_tax = d_gid_to_true_taxonomy[gid]
        true_tax_no_sp = true_tax[:true_tax.index('s__') - 1]

        # 1. Trivial case (everything but the species matches) = correct
        if decorated_tax.startswith(true_tax_no_sp):
            correct.add(gid)
            continue

        # Otherwise, non-trivial case, this could be a polyphyletic group,
        # or jumped at a higher rank. Need to check the reference tree.
        ref_highest_node, ref_taxon = get_highest_named_node_from_node(ref_tree.find_node_with_taxon_label(gid))
        dec_highest_node, dec_taxon = get_highest_named_node_from_node(dec_tree.find_node_with_taxon_label(raw_gid))

        # The named node should be the same, and leaf nodes should be the same (minus test gids)
        ref_desc_taxa = {x.taxon.label for x in ref_highest_node.leaf_nodes() if x.taxon.label not in test_gids}
        dec_desc_taxa = {x.taxon.label for x in dec_highest_node.leaf_nodes() if not x.taxon.label.startswith('TEST_')}

        if ref_taxon != dec_taxon or ref_desc_taxa != dec_desc_taxa or len(ref_desc_taxa) == 0:
            print('?')
            wrong.add(gid)
        else:
            correct.add(gid)

    return


def get_highest_named_node_from_node(node: dendropy.Node) -> Tuple[dendropy.Node, str]:

    parent_node = node.parent_node
    while parent_node is not None:

        # Parse the label information
        support, taxon, auxiliary_info = parse_label(parent_node.label)

        if taxon is not None and not taxon.startswith('s__'):
            return parent_node, taxon

        # Keep going up
        parent_node = parent_node.parent_node

    # Raise an exception if we can't find any named nodes
    raise Exception(f'Error finding named node for leaf: {node.label}')


def _read_decorated_tax_file(path):
    d_gid_to_decorated_tax = dict()
    with open(path) as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            d_gid_to_decorated_tax[gid] = tax.replace('; ', ';')
    return d_gid_to_decorated_tax


def parse_tax_string_into_dict_list(tax_string: str):
    """
    This function takes a taxonomic string and parses it into a list of
    dictionaries, where each dictionary represents a taxonomic rank.
    """
    d_rank_to_tax = defaultdict(list)
    for rank in tax_string.split(';'):
        d_rank_to_tax[rank[0]].append(rank)
    return d_rank_to_tax




def _test():
    from workflow.external.gtdb_metadata import GtdbMetadataR207
    from workflow.fasttree_full_tree_and_failed_no_trim.e_decorate_fasttree import FastTreeFullTreeAndFailedNoTrimDecorate
    from workflow.marker_trimming.gunc_e_decorate_fasttree import MarkerTrimmingGuncDecorateFastTree
    from workflow.util.log import log
    from workflow.fasttree_marker_analysis.gunc_d_analyse_decorated import FastTreeMarkerAnalyseDecorated


    log('Loading metadata')
    df_meta = GtdbMetadataR207().output().read_cached()

    log('Loading ref tree')
    ref_tree = FastTreeFullTreeAndFailedNoTrimDecorate().output().read()

    log('Loading dec tree')
    # dec_tree = MarkerTrimmingGuncDecorateFastTree(target_pct=50).output().read()
    dec_tree = FastTreeMarkerAnalysisDecorateGuncFastTree(target_pct=50, congruence=100).output().read()
    path_decorated_tax = '/srv/home/uqamussi/projects/gunc-chimeras/output/fasttree_marker_analysis/c100_p50/phylorank/decorated.tree-taxonomy'

    log('Comparing')
    check_phylorank_decorate_output(
        df_meta,
        ref_tree,
        dec_tree,
        path_decorated_tax
    )

    return


if __name__ == '__main__':
    _test()
