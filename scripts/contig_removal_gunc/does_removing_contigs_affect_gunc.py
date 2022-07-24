import os
import re
import shutil
from collections import defaultdict
from Bio import SeqIO
from magna.diamond import read_diamond_output
from magna.gunc import read_contig_assignments_tsv

REF_PATH = '/srv/home/uqamussi/projects/gunc-chimeras/luigi/genomes/GCA/000/255/595/prodigal/GCA_000255595.2.faa'
GENOME_DIR = '/srv/home/uqamussi/tmp/contig_removal_gunc/genomes'

diamond_root = '/srv/home/uqamussi/tmp/contig_removal_gunc/output/diamond_output'
gunc_root = '/srv/home/uqamussi/tmp/contig_removal_gunc/output/gunc_output'

import random

def read_fasta(path):
    out = dict()
    with open(path, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            out[str(record.description)] = str(record.seq)
    return out


def generate_data():
    os.makedirs(GENOME_DIR, exist_ok=True)

    faa = read_fasta(REF_PATH)

    contig_to_gene = defaultdict(set)
    contig_to_size = defaultdict(lambda: 0)
    for gene, seq in faa.items():
        contig = gene.split('_')[0]
        contig_to_gene[contig].add(gene)
        contig_to_size[contig] += len(seq)
    contig_to_size = sorted(contig_to_size.items(), key=lambda x: x[1], reverse=True)

    shutil.copy(REF_PATH, os.path.join(GENOME_DIR, 'GCA_000255595.2.faa'))

    faa_keys_in_order = tuple(faa.keys())
    for i in range(30):
        contigs_to_keep = set([x[0] for x in contig_to_size[i:]])
        path = os.path.join(GENOME_DIR, f'GCA_000255595.2_{i}_{contig_to_size[i][1]}.faa')
        tots_size = 0
        with open(path, 'w') as f:
            for gene, seq in faa.items():
                contig = gene.split('_')[0]
                if contig in contigs_to_keep:
                    f.write(f'>{gene}\n{seq}\n')
                    tots_size += len(seq)
        print(f'{path} {tots_size:,}')


def analyse_contig(gids):
    d_contig_to_cnt = dict()
    for gid in gids:
        contig_path = os.path.join(gunc_root, f'{gid}.contig_assignments.tsv')
        df_contig = read_contig_assignments_tsv(contig_path)
        for _, row in df_contig.iterrows():
            contig_key = (row['contig'], row['tax_level'], row['assignment'])
            new_val = row['count_of_genes_assigned']
            if contig_key not in d_contig_to_cnt:
                d_contig_to_cnt[contig_key] = new_val
            elif d_contig_to_cnt[contig_key] != new_val:
                print(f'{gid} {contig_key} {d_contig_to_cnt[contig_key]} != {new_val}')
            else:
                pass
    return


def analyse_diamond(gids):
    d_contig_to_cnt = dict()

    for gid in gids:
        diamond_path = os.path.join(diamond_root, f'{gid}.diamond.gtdb_95.out')

        df_diamond = read_diamond_output(diamond_path)

        for _, row in df_diamond.iterrows():
            k = tuple([row['query'], row['reference']])
            new_v = tuple(x for x in row.values[2:])
            if k not in d_contig_to_cnt:
                d_contig_to_cnt[k] = new_v
            elif d_contig_to_cnt[k] != new_v:
                print(f'{gid} {k} {d_contig_to_cnt[k]} != {new_v}')
            else:
                pass

    return


def analyse_data():
    gids = {x.replace('.diamond.gtdb_95.out', '') for x in os.listdir(diamond_root) if not x.startswith('._')}

    analyse_contig(gids)
    analyse_diamond(gids)

    return


def main():
    generate_data()
    # analyse_data()

    return


if __name__ == '__main__':
    main()
