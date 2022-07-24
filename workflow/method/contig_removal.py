import os
from collections import Counter
from collections import defaultdict
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from magna.gunc import read_contig_assignments_tsv

from workflow.config import PCT_VALUES
from workflow.util.fasta import read_fasta
from workflow.util.paths import get_gid_root


def get_taxonomy_by_majority_vote_gunc(df_contig_assign: pd.DataFrame, max_css: str, expected_domain: str) -> Tuple[
    str, str]:
    # Create a mapping from each tax_level to the contigs present
    d_contig_to_tax_level = defaultdict(lambda: defaultdict(dict))
    for _, row in df_contig_assign.iterrows():
        d_contig_to_tax_level[row['contig']][row['tax_level']][row['assignment']] = int(row['count_of_genes_assigned'])

    # Iterate over each taxonomic level until a majority has been reached
    majority_vote, majority_tax_level = None, None
    tax_levels = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus')
    for tax_level in tax_levels:

        # Take the maximum hit for each contig at this tax_level
        d_contig_to_vote = dict()
        for contig, d_tax_level_to_counts in d_contig_to_tax_level.items():
            cur_contig_assignments = d_tax_level_to_counts[tax_level]
            top_hit = sorted(cur_contig_assignments.items(), key=lambda x: x[1], reverse=True)[0]
            d_contig_to_vote[contig] = top_hit[0]

        # Aggregate the taxa counts
        d_taxon_to_count = Counter(d_contig_to_vote.values())
        total_count = sum(d_taxon_to_count.values())

        # Sanity check
        if tax_level == 'kingdom' and expected_domain not in d_taxon_to_count:
            raise Exception('This is a particularly bad case that you need to look into!')

        # Do a majority vote
        majority_rank, majority_count = d_taxon_to_count.most_common(1)[0]
        if majority_count / total_count > 0.5:
            majority_tax_level = tax_level
            if tax_level == 'kingdom' and majority_rank != expected_domain:
                majority_vote = expected_domain
                break
            else:
                majority_vote = majority_rank
        else:
            # Unable to get a majority at the domain level, force the domain to be the expected domain and stop
            if tax_level == 'kingdom':
                majority_vote = expected_domain
                majority_tax_level = tax_level
                break

            # No majority vote exists, just use the previous iteration majority values
            else:
                break

        # Stop processing if we have reached the MaxCSS level reported by GUNC
        if tax_level == max_css:
            break

    if majority_vote is None or majority_tax_level is None:
        raise Exception('No majority vote found!')

    return majority_vote, majority_tax_level


def get_contig_metadata(d_fna: Dict[str, str]):
    out = dict()
    for contig, seq in d_fna.items():
        seq = str(seq.seq)
        contig_len = len(seq)
        gc_pct = 100 * (seq.count('G') + seq.count('C')) / contig_len
        out[contig] = {'length': len(seq), 'gc': gc_pct}
    return out


def contigs_to_remove_from_gunc(d_fna: Dict[str, str], df_contig_assign: pd.DataFrame, taxon: str, tax_level: str):
    # Subset the dataframe to only those ranks in the tax level
    df_subset = df_contig_assign[df_contig_assign['tax_level'] == tax_level]
    if len(df_subset) == 0:
        raise Exception(f'No contigs found at tax_level {tax_level}')

    # Aggregate the assignments by contig
    d_contig_to_assignments = defaultdict(dict)
    for _, row in df_subset.iterrows():
        d_contig_to_assignments[row['contig']][row['assignment']] = int(row['count_of_genes_assigned'])

    # Get the mapping percentage for each contig
    d_contig_to_map_pct = dict()
    for contig, d_assignments in d_contig_to_assignments.items():
        total_count = sum(d_assignments.values())
        d_contig_to_map_pct[contig] = 100 * d_assignments.get(taxon, 0) / total_count

    # Load the contig metadata
    d_contig_metadata = get_contig_metadata(d_fna)

    # Contigs which have >50% of the majority vote can be considered for the GC% calculation
    contigs_for_gc_calculation = {c for c, p in d_contig_to_map_pct.items() if p > 50}
    d_contig_to_gc_delta = dict()
    if len(contigs_for_gc_calculation) > 0:
        median_gc = float(np.median([d_contig_metadata[x]['gc'] for x in contigs_for_gc_calculation]))
        for contig, contig_meta in d_contig_metadata.items():
            d_contig_to_gc_delta[contig] = abs(contig_meta['gc'] - median_gc)

    # Create an ordering from most to least likely to be contaminated
    # Sort by %map (asc), length (asc), delta GC% (desc)
    contigs_to_sort = list()
    for contig, map_pct in d_contig_to_map_pct.items():
        contigs_to_sort.append((
            contig,
            map_pct,
            d_contig_metadata[contig]['length'],
            d_contig_to_gc_delta.get(contig, 0)
        ))
    contig_ordering = sorted(contigs_to_sort, key=lambda x: (x[1], x[2], -x[3]))

    # Create the percentages at which contigs can be removed
    genome_len = sum([d_contig_metadata[x]['length'] for x in d_contig_metadata])
    d_pct_to_contigs = dict()
    for pct in PCT_VALUES:
        contigs_to_remove = set()
        total_len_removed = 0
        for contig, _, contig_len, _ in contig_ordering:
            new_len_removed = total_len_removed + contig_len
            new_pct_removed = new_len_removed / genome_len * 100

            # Don't add it as it will be over the threshold
            if new_pct_removed > pct:
                continue

            contigs_to_remove.add(contig)
            total_len_removed += contig_len

        d_pct_to_contigs[pct] = contigs_to_remove

    # De-duplicate the data set to only those percentages where a unique number of contigs were removed
    last_contigs_removed = set()
    d_pct_to_contigs_dedup = dict()
    for pct, contigs_removed in d_pct_to_contigs.items():
        if contigs_removed == last_contigs_removed:
            continue
        d_pct_to_contigs_dedup[pct] = contigs_removed
        last_contigs_removed = contigs_removed
    return d_pct_to_contigs_dedup


def _main():
    gid = 'GCF_003052605.1'
    domain = 'd__Bacteria'
    max_css = 'genus'

    gid_root = get_gid_root(gid)

    # Read the fasta file
    d_fna = read_fasta(os.path.join(gid_root, f'{gid}.fna'))

    gunc_dir_r95 = os.path.join(gid_root, 'gunc_r95')
    gunc_dir_pro = os.path.join(gid_root, 'gunc_pro')

    df_contig_assign = read_contig_assignments_tsv(
        os.path.join(gunc_dir_r95, 'gunc_output', f'{gid}.contig_assignments.tsv'))

    taxon, tax_level = get_taxonomy_by_majority_vote_gunc(df_contig_assign, max_css, domain)
    d_contigs_to_remove = contigs_to_remove_from_gunc(d_fna, df_contig_assign, taxon, tax_level)

    return


if __name__ == '__main__':
    _main()
