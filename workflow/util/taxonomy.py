from collections import defaultdict
from typing import Dict, Set, Tuple

D_TAX_NOVELTY_SUFFIX_TO_NAME = {
    'd': 'domain',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
    'st': 'strain'
}

D_TAX_NOVELTY_NAME_TO_SUFFIX = {
    'domain': 'd',
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's',
    'strain': 'st'
}

D_TAX_NOVELTY_TO_INDEX = {
    'domain': 0,
    'phylum': 1,
    'class': 2,
    'order': 3,
    'family': 4,
    'genus': 5,
    'species': 6,
    'strain': 7
}

INDEX_TO_TAX_NOVELTY = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')


def calculate_taxonomic_novelty(d_gid_to_tax_string: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    """Calculates the degree of taxonomic novelty for a collection of taxa."""
    d_tax_under = defaultdict(lambda: 0)
    for taxonomy in d_gid_to_tax_string.values():
        for rank in taxonomy.replace('; ', ';').split(';'):
            d_tax_under[rank] += 1

    # Find the novelty
    d_gid_to_tax_novelty = dict()
    d_tax_novelty_to_gid = dict()
    for gid, tax_str in d_gid_to_tax_string.items():
        taxonomy = tax_str.replace('; ', ';').split(';')
        last_rank = 'st'
        for taxon in reversed(taxonomy):
            taxon_count = d_tax_under[taxon]
            if taxon_count > 1:
                break
            last_rank = taxon[0]

        tax_novelty = D_TAX_NOVELTY_SUFFIX_TO_NAME[last_rank]

        d_gid_to_tax_novelty[gid] = tax_novelty

        if tax_novelty not in d_tax_novelty_to_gid:
            d_tax_novelty_to_gid[tax_novelty] = set()
        d_tax_novelty_to_gid[tax_novelty].add(gid)

    return d_gid_to_tax_novelty, d_tax_novelty_to_gid
