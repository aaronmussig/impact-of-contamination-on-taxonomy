
import os
from collections import defaultdict

from tqdm import tqdm

from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.method.run_mummer_on_genome import COORD_HEADER
from workflow.util.fasta import read_fasta
from workflow.util.paths import get_gid_root
import pandas as pd
import numpy as np


GID = 'GCA_000153745.1'
PATH = os.path.join(get_gid_root(GID), 'nucmer_to_rep_genomes.h5')
PATH_FNA = os.path.join(get_gid_root(GID), f'{GID}.fna')
# PATH = '/srv/home/uqamussi/projects/gunc-chimeras/output/genomes/GCA/002/170/165/nucmer_to_rep_genomes.h5'

def parse_for_gid(q_gid, rows):
    df = pd.DataFrame(rows, columns=COORD_HEADER)

    est_len = 0
    for row in df.itertuples():
        est_len += row.len_1
    return est_len

def get_fna_contigs_and_lengths(path_fna):
    out = dict()
    for contig, seq in read_fasta(path_fna).items():
        out[contig] = len(seq)
    return out

def group_by_contig_to_gid_rows(df):
    out_tmp = defaultdict(list)
    for row in df.itertuples():
        out_tmp[row.contig_r].append(row[1:])

    out = dict()
    for contig_r, rows in out_tmp.items():
        df_new = pd.DataFrame(rows, columns=COORD_HEADER)
        out[contig_r] = df_new
    return out


def get_coverage_for_qry_gids(df, contig_len):
    d_gid_to_hits = defaultdict(lambda: 0)
    for row in tqdm(df.itertuples(), total=len(df)):
        d_gid_to_hits[row.qry_gid] += row.len_1 * row.pct_identity
    return d_gid_to_hits



def parse_file(path):

    d_contig_to_len = get_fna_contigs_and_lengths(PATH_FNA)

    df_meta = GtdbMetadataR207().output().read_cached()

    r_tax = df_meta.loc[GID, 'gtdb_taxonomy']
    print(f'{GID}\t{r_tax}')

    df = pd.read_hdf(path)

    d_contig_to_df = group_by_contig_to_gid_rows(df)

    for contig_r, contig_df in d_contig_to_df.items():
        print(d_contig_to_len[contig_r])

        qry_gid_to_cov = get_coverage_for_qry_gids(contig_df, d_contig_to_len[contig_r])

        i = 0
        for qry_gid, cov in sorted(qry_gid_to_cov.items(), key=lambda x: x[1], reverse=True):
            print(qry_gid, cov, df_meta.loc[qry_gid, 'gtdb_taxonomy'])
            i += 1
            if i > 100:
                break

    return


    d_gid_to_hits = defaultdict(list)
    for row in tqdm(df.itertuples(), total=len(df), smoothing=0.01):
        d_gid_to_hits[row.qry_gid].append(tuple(row[1:]))

    q_gid_to_cnt = dict()
    for q_gid, rows in tqdm(d_gid_to_hits.items(), smoothing=0.01):
        gid_results = parse_for_gid(q_gid, rows)
        q_gid_to_cnt[q_gid] = gid_results

    i = 0

    for q_gid, cnt in sorted(q_gid_to_cnt.items(), key=lambda x: x[1], reverse=True):
        i += 1
        q_tax = df_meta.loc[q_gid, 'gtdb_taxonomy']
        print(f'{q_gid}\t{cnt}\t{q_tax}')
        if i > 100:
            break

    return




if __name__ == '__main__':
    parse_file(PATH)