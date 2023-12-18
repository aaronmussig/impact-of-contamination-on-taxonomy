import multiprocessing as mp
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_GUNC
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.get_taxonomy_from_gunc_contig_assignment import GetTaxonomyFromGuncContigAssignment
from workflow.method.contig_removal import get_contig_metadata
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class GuncRankContigsByContamination(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'contig_assignments': GetTaxonomyFromGuncContigAssignment()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_GUNC, 'gunc_contig_ranking.h5'))

    def run(self):
        log('Loading contig assignments')
        df_contig_assignments = self.input()['contig_assignments'].maybe_read_cached()

        log('Loading max CSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()

        log('Creating queue')
        queue = list()
        for gid, row in tqdm(df_max_css.iterrows(), total=len(df_max_css)):
            inferred_tax = df_contig_assignments.loc[gid, 'taxonomy']
            max_css = row['taxonomic_level']
            queue.append((gid, inferred_tax, max_css, row['source']))

            if DEBUG and len(queue) > 10:
                break


        with mp.Pool(processes=mp.cpu_count() ) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        df = pd.concat(results)
        df.sort_values(by=['gid', 'order'], inplace=True, ignore_index=True)
        #
        if not DEBUG:
            self.save_hdf(df)
        return


RANKS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


def worker(job):
    gid, tax, max_css, source = job

    # Calculate the highest achieved clade separation score when inferring taxonomy
    inferred_max_css = RANKS[len(tax.split(';')) - 1]
    if RANKS.index(max_css) > RANKS.index(inferred_max_css):
        raise Exception(f'{gid} has max CSS {max_css} but inferred taxonomy has {inferred_max_css}')

    gid_root = get_gid_root(gid)
    d_fna = read_fasta(os.path.join(gid_root, f'{gid}.fna'))

    if source == 'progenomes':
        gunc_dir = os.path.join(gid_root, 'gunc_pro')
    elif source == 'gtdb':
        gunc_dir = os.path.join(gid_root, 'gunc_r95')
    else:
        raise ValueError(f'Unknown source: {source}')

    df_contig_assign = read_contig_assignments_tsv(
        os.path.join(gunc_dir, 'gunc_output', f'{gid}.contig_assignments.tsv'))
    df_subset = df_contig_assign[df_contig_assign['tax_level'] == max_css]
    taxon = tax.split(';')[RANKS.index(max_css)]

    # Calculate the percentage of genes on each contig that mapped to the taxon
    d_contig_to_total_counts = defaultdict(lambda: 0)
    d_contig_to_correct_counts = defaultdict(lambda: 0)
    for _, row in df_subset.iterrows():
        if row['assignment'] == taxon:
            d_contig_to_correct_counts[row['contig']] += int(row['count_of_genes_assigned'])
        d_contig_to_total_counts[row['contig']] += int(row['count_of_genes_assigned'])

    d_contig_metadata = get_contig_metadata(d_fna)

    out = list()
    for contig, contig_metadata in d_contig_metadata.items():
        total_gene_cnt = d_contig_to_total_counts.get(contig, 0)
        correct_gene_cnt = d_contig_to_correct_counts.get(contig, 0)
        pct_correct = 100 * (correct_gene_cnt / total_gene_cnt) if total_gene_cnt > 0 else 0
        out.append({
            'gid': gid,
            'contig': contig,
            'length': contig_metadata['length'],
            'gc': contig_metadata['gc'],
            'correct_genes': correct_gene_cnt,
            'total_genes': total_gene_cnt,
            'pct_correct': pct_correct
        })

    df = pd.DataFrame(out)

    # Get the GC% of the contig that has the largest number of genes correct
    row_largest_correct = df.sort_values(by='correct_genes', ascending=False).iloc[0]
    gc_largest_correct = row_largest_correct['gc']
    df['gc_delta'] = abs(df['gc'] - gc_largest_correct)

    # Create the ordering
    df.sort_values(by=['pct_correct', 'length', 'gc_delta'], ascending=[True, True, False], inplace=True)
    df['order'] = np.arange(0, len(df))
    return df
