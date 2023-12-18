from collections import deque, defaultdict
from typing import Dict, Set
from typing import FrozenSet

import dendropy
from tqdm import tqdm

from workflow.util.obj import is_float


def parse_label(label):
    """Parse a Newick label which may contain a support value, taxon, and/or auxiliary information.

    Parameters
    ----------
    label : str
        Internal label in a Newick tree.

    Returns
    -------
    float
        Support value specified by label, or None
    str
        Taxon specified by label, or None
    str
        Auxiliary information, on None
    """

    support = None
    taxon = None
    auxiliary_info = None

    if label:
        label = label.strip()
        if '|' in label:
            label, auxiliary_info = label.split('|')

        if ':' in label:
            support, taxon = label.split(':')
            support = float(support)
        else:
            if is_float(label):
                support = float(label)
            elif label != '':
                taxon = label

    return support, taxon, auxiliary_info


def get_tree_node_depth(tree: dendropy.Tree) -> Dict[int, Set[dendropy.Node]]:
    """Calculate the depth of all nodes in the tree.

    :param tree: The dendropy tree object to calculate the depths for.

    :returns: A dictionary of depths to nodes at that depth.
    """
    # Validate arguments
    if not isinstance(tree, dendropy.Tree):
        raise ValueError(f'Tree is not a dendropy.Tree object.')

    # Create a queue starting from the seed node (depth=0)
    out = defaultdict(set)
    queue = deque([(tree.seed_node, 0)])

    # Progressively add children to the queue
    while len(queue) > 0:
        cur_node, cur_depth = queue.popleft()
        out[cur_depth].add(cur_node)

        for child in cur_node.child_nodes():
            queue.append((child, cur_depth + 1))

    # Make this a dictionary (not defaultdict)
    return dict(out)


def get_tree_node_to_desc_taxa(tree: dendropy.Tree) -> Dict[dendropy.Node, FrozenSet[str]]:
    """Get each node and the taxa that are descendants of this node.

    :param tree: The dendropy tree object to calculate the depths for.

    :returns: A dictionary of nodes and their descendant taxa.
    """

    # Calculate the nodes at each depth
    node_depth = get_tree_node_depth(tree)

    # Process the deepest nodes first
    out = defaultdict(set)
    for depth, nodes in sorted(node_depth.items(), key=lambda x: x[0], reverse=True):
        for node in nodes:

            # Seed the output dictionary
            if node.is_leaf():
                out[node].add(node.taxon.label)

            # Check bring up the taxa from the descendant nodes
            else:
                for child in node.child_nodes():
                    out[node].update(out[child])

    assert (out[tree.seed_node] == {x.label for x in tree.taxon_namespace})
    return {k: frozenset(v) for k, v in out.items()}


def read_taxonomy_from_tree(tree: dendropy.Tree) -> Dict[str, str]:
    d_node_to_tax_strings = dict()

    for node in tqdm(tree.leaf_node_iter(), total=len(tree.taxon_namespace)):
        tax_strings = defaultdict(list)

        cur_node = node
        while cur_node is not None:
            support, taxa, _ = parse_label(cur_node.label)
            if taxa is not None:
                for taxon in taxa.replace('; ', ';').split(';'):
                    tax_strings[taxon[0]].append(taxon)
            cur_node = cur_node.parent_node

        d_node_to_tax_strings[node.taxon.label] = tax_strings

    # Convert into a proper tax string, and fill in blank ranks
    out = dict()
    for gid, d_tax_strings in d_node_to_tax_strings.items():
        cur_tax_string = list()
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            if rank in d_tax_strings:
                for taxon in d_tax_strings[rank]:
                    cur_tax_string.append(taxon)
            else:
                cur_tax_string.append(f'{rank}__')
        out[gid] = ';'.join(cur_tax_string)
    return out


def read_taxonomy_from_tree_pplacer(tree: dendropy.Tree, d_gid_to_tax) -> Dict[str, str]:
    d_node_to_tax_strings = dict()

    for node in tree.leaf_node_iter():

        # Only interested in test genomes
        if not node.taxon.label.startswith('TEST_'):
            continue

        # If this node a sister to a reference genome, then take the taxonomy from the reference genome
        v_cur_node = node
        skip_node = False
        while v_cur_node is not None:
            sister_nodes = v_cur_node.leaf_nodes()
            non_test_gids = [not x.taxon.label.startswith('TEST_') for x in sister_nodes]

            # Try get a tax string
            if any(non_test_gids):
                sister_taxonomy_strings = set()
                for sister_node in sister_nodes:
                    if not sister_node.taxon.label.startswith('TEST_'):
                        sister_tax_string = d_gid_to_tax[sister_node.taxon.label[3:]]
                        sister_taxonomy_strings.add(';'.join(sister_tax_string.split(';')[0:6]))
                if len(sister_taxonomy_strings) == 1:
                    d_node_to_tax_strings[node.taxon.label] = {y[0]: [y] for y in  sister_taxonomy_strings.pop().split(';')}
                    skip_node = True

                # Break from the while loop, since there is a disagreement, must read from the tree
                break

            # otherwise, keep going up.
            v_cur_node = v_cur_node.parent_node

        # We found the tax using a sister, don't search the tree
        if skip_node:
            continue

        sister_nodes = node.parent_node.leaf_nodes()
        sister_taxonomy_strings = set()
        for sister_node in sister_nodes:
            if not sister_node.taxon.label.startswith('TEST_'):
                sister_tax_string = d_gid_to_tax[sister_node.taxon.label[3:]]
                sister_taxonomy_strings.add(';'.join(sister_tax_string.split(';')[0:6]))
        if len(sister_taxonomy_strings) == 1:
            d_node_to_tax_strings[node.taxon.label] = {y[0]: [y] for y in sister_taxonomy_strings.pop().split(';')}
            continue
        else:
            print()

        # Otherwise, infer from the tree
        tax_strings = defaultdict(list)

        cur_node = node
        while cur_node is not None:
            support, taxa, _ = parse_label(cur_node.label)
            if taxa is not None:
                for taxon in taxa.replace('; ', ';').split(';'):
                    tax_strings[taxon[0]].append(taxon)
            cur_node = cur_node.parent_node

        d_node_to_tax_strings[node.taxon.label] = tax_strings

    # Convert into a proper tax string, and fill in blank ranks
    out = dict()
    for gid, d_tax_strings in d_node_to_tax_strings.items():
        cur_tax_string = list()
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            if rank in d_tax_strings:
                for taxon in d_tax_strings[rank]:
                    cur_tax_string.append(taxon)
            else:
                cur_tax_string.append(f'{rank}__')
        out[gid] = ';'.join(cur_tax_string)
    return out