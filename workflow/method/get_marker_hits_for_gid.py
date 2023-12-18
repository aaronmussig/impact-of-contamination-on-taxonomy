import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Set

from workflow.config import R207_MARKERS
from workflow.method.get_genome_domain_from_markers import get_genome_domain_from_markers
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile, Hit
from workflow.util.fasta import read_fasta
from workflow.util.paths import get_gid_root


def _merge_hit_files(pfam_th: TopHitPfamFile, tigr_th: TopHitTigrFile, omit_contigs: Set[str]) -> Dict[str, List[Hit]]:
    """Returns a list of Hits grouped by HMM id as a dictionary."""
    out = dict()
    for cur_tophit in (pfam_th, tigr_th):
        for gene_id, cur_hit in cur_tophit.iter_hits():
            contig_id = gene_id[:gene_id.rindex('_')]
            if contig_id not in omit_contigs:
                if cur_hit.hmm_id not in out:
                    out[cur_hit.hmm_id] = list()
                out[cur_hit.hmm_id].append(cur_hit)
    return out



def get_all_marker_hits_for_gid(d_genes, pfam_th, tigr_th):
    """Return all and their seq"""

    # Pointers to unique, multiple hit, multiple-unique, missing markers.
    cur_unq = dict()
    cur_mul = dict()
    cur_muq = dict()
    cur_mis = dict()

    re_pos = re.compile(r'.+? # (\d+) # (\d+) #.+')

    d_marker_to_position = defaultdict(list)

    # Load genes from the prodigal faa file.
    for seq_id, seq in d_genes.items():
        if str(seq.seq).endswith('*'):
            d_genes[seq_id].seq = seq.seq[:-1]

    # Create a dictionary of marker names -> Hits
    d_hmm_hits = _merge_hit_files(pfam_th, tigr_th, set())

    # Foreach expected marker determine which category it falls into.
    for marker_id in R207_MARKERS:

        # Marker is missing.
        if marker_id not in d_hmm_hits:
            cur_mis[marker_id] = None

        # Multiple hits to the same marker.
        elif len(d_hmm_hits[marker_id]) > 1:

            # If sequences are the same, take the most significant hit
            unq_seqs = {str(d_genes[x.gene_id].seq) for x in d_hmm_hits[marker_id]}
            if len(unq_seqs) == 1:
                cur_top_hit = sorted(d_hmm_hits[marker_id], reverse=True)[0]
                cur_muq[marker_id] = {'hit': cur_top_hit, 'seq': d_genes[cur_top_hit.gene_id]}
                print("NOT IMPLEMENTED!")
                for marker_hit in d_hmm_hits[marker_id]:
                    cur_hit = d_hmm_hits[marker_id][0]
                    cur_seq = d_genes[cur_hit.gene_id]

                    re_hit = re_pos.search(cur_seq.description)
                    if re_hit is None:
                        raise Exception('??')
                    pos_from, pos_to = re_hit.groups()
                    pos_from, pos_to = int(pos_from), int(pos_to)
                    d_marker_to_position[marker_id].append((cur_hit.gene_id, pos_from, pos_to))


            # Marker maps to multiple genes.
            else:
                for cur_hit in d_hmm_hits[marker_id]:
                    cur_seq = d_genes[cur_hit.gene_id]

                    re_hit = re_pos.search(cur_seq.description)
                    if re_hit is None:
                        raise Exception('??')
                    pos_from, pos_to = re_hit.groups()
                    pos_from, pos_to = int(pos_from), int(pos_to)
                    d_marker_to_position[marker_id].append((cur_hit.gene_id, pos_from, pos_to))

                cur_mul[marker_id] = None

        # This was a unique hit.
        else:
            cur_hit = d_hmm_hits[marker_id][0]
            cur_seq = d_genes[cur_hit.gene_id]

            re_hit = re_pos.search(cur_seq.description)
            if re_hit is None:
                raise Exception('??')
            pos_from, pos_to = re_hit.groups()
            pos_from, pos_to = int(pos_from), int(pos_to)
            d_marker_to_position[marker_id].append((cur_hit.gene_id, pos_from, pos_to))

            cur_unq[marker_id] = {'hit': cur_hit, 'seq': d_genes[cur_hit.gene_id]}

    out_markers = {
        'unq': cur_unq,
        'mul': cur_mul,
        'muq': cur_muq,
        'mis': cur_mis,
    }
    return d_marker_to_position, out_markers

def get_marker_hits_for_gid(d_genes, pfam_th, tigr_th, omit_contigs: Optional[Set[str]] = None):
    if omit_contigs is None:
        omit_contigs = set()

    # Pointers to unique, multiple hit, multiple-unique, missing markers.
    cur_unq = dict()
    cur_mul = dict()
    cur_muq = dict()
    cur_mis = dict()

    # Load genes from the prodigal faa file.
    for seq_id, seq in d_genes.items():
        if seq.endswith('*'):
            d_genes[seq_id] = seq[:-1]

    # Create a dictionary of marker names -> Hits
    d_hmm_hits = _merge_hit_files(pfam_th, tigr_th, omit_contigs)

    # Foreach expected marker determine which category it falls into.
    for marker_id in R207_MARKERS:

        # Marker is missing.
        if marker_id not in d_hmm_hits:
            cur_mis[marker_id] = None

        # Multiple hits to the same marker.
        elif len(d_hmm_hits[marker_id]) > 1:

            # If sequences are the same, take the most significant hit
            unq_seqs = {d_genes[x.gene_id] for x in d_hmm_hits[marker_id]}
            if len(unq_seqs) == 1:
                cur_top_hit = sorted(d_hmm_hits[marker_id], reverse=True)[0]
                cur_muq[marker_id] = {'hit': cur_top_hit, 'seq': d_genes[cur_top_hit.gene_id]}

            # Marker maps to multiple genes.
            else:
                cur_mul[marker_id] = None

        # This was a unique hit.
        else:
            cur_hit = d_hmm_hits[marker_id][0]
            cur_unq[marker_id] = {'hit': cur_hit, 'seq': d_genes[cur_hit.gene_id]}

    # invert the dict
    out_markers = {
        'unq': cur_unq,
        'mul': cur_mul,
        'muq': cur_muq,
        'mis': cur_mis,
    }
    return out_markers


def get_marker_hits_for_gid_include_mul(d_genes, pfam_th, tigr_th):
    omit_contigs = set()

    # Pointers to unique, multiple hit, multiple-unique, missing markers.
    cur_unq = dict()
    cur_mul = dict()
    cur_muq = dict()
    cur_mis = dict()

    out = list()

    # Load genes from the prodigal faa file.
    for seq_id, seq in d_genes.items():
        if seq.endswith('*'):
            d_genes[seq_id] = seq[:-1]

    # Create a dictionary of marker names -> Hits
    d_hmm_hits = _merge_hit_files(pfam_th, tigr_th, omit_contigs)

    # Foreach expected marker determine which category it falls into.
    for marker_id in R207_MARKERS:

        # Marker is missing.
        if marker_id not in d_hmm_hits:
            cur_mis[marker_id] = None

        # Multiple hits to the same marker.
        elif len(d_hmm_hits[marker_id]) > 1:

            # If sequences are the same, take the most significant hit
            unq_seqs = {d_genes[x.gene_id] for x in d_hmm_hits[marker_id]}
            if len(unq_seqs) == 1:
                cur_top_hit = sorted(d_hmm_hits[marker_id], reverse=True)[0]
                cur_muq[marker_id] = {'hit': cur_top_hit, 'seq': d_genes[cur_top_hit.gene_id]}
                out.append({'hit': cur_top_hit, 'n': 0, 'seq': d_genes[cur_top_hit.gene_id]})

            # Marker maps to multiple genes.
            else:
                cur_mul[marker_id] = None
                for x_i, x in enumerate(d_hmm_hits[marker_id]):
                    out.append({'hit': x, 'n': x_i, 'seq': d_genes[x.gene_id]})

        # This was a unique hit.
        else:
            cur_hit = d_hmm_hits[marker_id][0]
            cur_unq[marker_id] = {'hit': cur_hit, 'seq': d_genes[cur_hit.gene_id]}
            out.append({'hit': cur_hit, 'n': 0, 'seq': d_genes[cur_hit.gene_id]})

    # invert the dict
    out_markers = {
        'unq': cur_unq,
        'mul': cur_mul,
        'muq': cur_muq,
        'mis': cur_mis,
    }
    return out_markers, out


def main():
    gid = 'GCF_003052605.1'
    gid_root = get_gid_root(gid)

    # Read the fasta file
    d_faa = dict(read_fasta(os.path.join(gid_root, 'prodigal', f'{gid}.faa')))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    omit_contigs = {'NZ_QASO01000049.1'}
    # omit_contigs = None
    results_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th, omit_contigs)

    domain = get_genome_domain_from_markers(results_markers)

    return


if __name__ == '__main__':
    main()
