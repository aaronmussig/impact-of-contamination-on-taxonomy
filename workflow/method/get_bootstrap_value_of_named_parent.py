import dendropy
from tqdm import tqdm

from workflow.external.gtdb_tree import GtdbTreeBacR207
from workflow.util.tree import parse_label


def check_ref_tree(ref_tree: dendropy.Tree, gid_set, set_name, batch_id):
    out = list()

    for gid in tqdm(sorted(gid_set)):
        highest_taxon = None
        highest_bs = None

        leaf = ref_tree.find_node_with_taxon_label(gid)
        cur_parent = leaf.parent_node
        while cur_parent is not None:

            support, taxon, auxiliary_info = parse_label(cur_parent.label)

            if taxon is not None:
                highest_bs = support
                highest_taxon = taxon
                break

            cur_parent = cur_parent.parent_node

        if highest_bs is None or highest_taxon is None:
            out.append({
                'gid': gid,
                'batch_id': batch_id,
                'bs': -100,
                'set': set_name,
                'taxon': 'N/A',
            })

        else:
            out.append({
                'batch_id': batch_id,
                'set': set_name,
                'gid': gid,
                'bs': int(highest_bs),
                'taxon': highest_taxon,
            })


    return out



def main():
    ref_tree = GtdbTreeBacR207().output().read()
    for leaf_node in ref_tree.leaf_node_iter():
        leaf_node.taxon.label = leaf_node.taxon.label[3:]


    gid_set = {'GCA_002402315.1', 'GCA_003696645.1', 'GCA_012964215.1', 'GCA_013822275.1', 'GCA_018822835.1', 'GCA_019241735.1', 'GCF_000092025.1', 'GCF_002814095.1', 'GCF_002833405.1', 'GCF_003597675.1', 'GCF_010669145.1', 'GCF_013361095.1'}
    results = check_ref_tree(ref_tree, gid_set, 'asd', '0')
    return


if __name__ == '__main__':
    main()