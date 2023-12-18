import os
from random import shuffle

from Bio import SeqIO

from workflow.util.paths import get_gid_root


def randomly_select_contigs_up_to_pct(d_fna, pct):
    all_contigs = frozenset(d_fna.keys())

    # Get the length of all contigs
    d_contig_to_len = dict()
    for contig_id, seq in d_fna.items():
        d_contig_to_len[contig_id] = len(seq)
    genome_len = sum(d_contig_to_len.values())

    # Determine a random ordering
    contig_ordering = list(d_contig_to_len.keys())
    shuffle(contig_ordering)

    # Randomly select contigs to remove up to pct
    contigs_to_remove = set()
    total_len_removed = 0
    for contig in contig_ordering:
        contig_len = d_contig_to_len[contig]
        new_len_removed = total_len_removed + contig_len
        new_pct_removed = new_len_removed / genome_len * 100

        # Don't add it as it will be over the threshold
        if new_pct_removed > pct:
            continue

        contigs_to_remove.add(contig)
        total_len_removed += contig_len

    # Return those contigs to remove, and those to keep
    contigs_to_keep = all_contigs - contigs_to_remove
    return contigs_to_keep, contigs_to_remove


def main():
    gid = 'GCA_002170165.1'
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')

    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    randomly_select_contigs_up_to_pct(d_fna, 20)

    return


if __name__ == '__main__':
    main()
