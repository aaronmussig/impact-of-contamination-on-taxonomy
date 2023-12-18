from collections import Counter

import pandas as pd

from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class BlastGetTaxStringAgainstHits(LuigiTask):

    def requires(self):
        return {
            # 'blast': GtdbR207BlastDb(),
            # 'css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207()
        }

    # def output(self):
    #     return LocalTargetHdf5(os.path.join(DIR_OUT_BLAST, f'blast_genome_against_207_db.h5'))

    def run(self):
        log(f'BlastGetTaxStringAgainstHits', title=True)

        path = '/srv/home/uqamussi/tmp/blast_results.out'

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()

        df = parse_blast_output(path)

        top_hits = list()

        for q_gid, q_contig, df_subset in iter_dataframe_into_groups(df):

            q_tax = d_gid_to_tax[q_gid]
            # r_tax_counts = Counter([d_gid_to_tax[x] for x in df_subset['r_gid'].unique()])
            r_tax_counts = Counter([d_gid_to_tax[x] for x in df_subset[df_subset['pct_identity'] >= 90]['r_gid'].unique()])

            lst_counts = sorted(r_tax_counts.items(), key=lambda x: x[1], reverse=True)

            top_hits.append(lst_counts[0][0])

        total_counts = Counter(top_hits)
        for tax, count in total_counts.items():
            print(tax, count)
        return


def iter_dataframe_into_groups(df):
    # sort the dataframe
    df.sort_values(
        by=['q_gid', 'q_contig', 'bit_score', 'pct_identity', 'aln_len'],
        ascending=[True, True, False, False, False],
        inplace=True,
        ignore_index=True
    )

    last_key = None
    last_i = 0
    for i in range(len(df)):
        cur_row = df.iloc[i]
        cur_key = (cur_row.q_gid, cur_row.q_contig)
        if last_key is None:
            last_key = cur_key
        elif last_key != cur_key:

            # Subset the dataframe and yield the results
            df_subset = df.iloc[last_i:i]
            yield last_key[0], last_key[1], df_subset

            # Store the new value and index
            last_key = cur_key
            last_i = i

    return


def read_genome_fna(gid):
    genome = Genome(gid)
    d_fna = genome.d_fna
    return gid, d_fna


def parse_blast_output(path):
    rows = list()

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')

            q_gid, q_contig = cols[0].split('__')
            s_gid, s_contig = cols[1].split('__')

            rows.append({
                'q_gid': q_gid,
                'q_contig': q_contig,
                'r_gid': s_gid,
                'r_contig': s_contig,
                'pct_identity': float(cols[2]),
                'aln_len': int(cols[3]),
                'mismatches': int(cols[4]),
                'gap_openings': int(cols[5]),
                'q_start': int(cols[6]),
                'q_end': int(cols[7]),
                'r_start': int(cols[8]),
                'r_end': int(cols[9]),
                'e_value': float(cols[10]),
                'bit_score': float(cols[11]),
            })

    df = pd.DataFrame(rows)

    return df
