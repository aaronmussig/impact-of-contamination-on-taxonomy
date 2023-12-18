import os
from collections import defaultdict
import multiprocessing as mp

import pandas as pd
from gunc.get_scores import read_diamond_output, create_base_data, get_stats, get_abundant_lineages_cutoff, \
    get_scores_for_taxlevel
from magna.gunc import read_contig_assignments_tsv
from itertools import combinations

from tqdm import tqdm

from workflow.config import TAX_LEVELS
from workflow.external.gtdb_metadata import GtdbMetadataR95
from workflow.gunc_helper.aggregate_max_css_level_gunc import AggregateMaxCssLevelGtdbR95
from workflow.util.fasta import read_fasta
from workflow.util.gunc import get_contig_mapping_proportion_for_taxon
from workflow.util.paths import get_gid_root


def get_df():
    df_meta = GtdbMetadataR95().output().read_cached()
    # df_max_css = AggregateMaxCssLevelProGenomes().output().read_cached()
    df_max_css = AggregateMaxCssLevelGtdbR95().output().read_cached()
    df = pd.merge(df_meta, df_max_css, left_index=True, right_index=True)
    df = df[df['pass.GUNC'] == False]
    return df


def get_taxonomy_by_majority(path_contig_assign, taxonomic_level):
    df_contig_assign = read_contig_assignments_tsv(path_contig_assign)
    df_contig_assign = df_contig_assign[df_contig_assign['tax_level'] == taxonomic_level]

    contig_to_tax = defaultdict(list)
    for _, row in df_contig_assign.iterrows():
        contig_to_tax[row['contig']].append((row['assignment'], row['count_of_genes_assigned']))

    d_tax_to_count = defaultdict(lambda: 0)
    for contig, lst_tax_assignments in contig_to_tax.items():
        d_tax_to_cnt = dict()
        for tax, cnt in lst_tax_assignments:
            d_tax_to_cnt[tax] = cnt
            d_tax_to_count[tax] += cnt

    winner = sorted(d_tax_to_count.items(), key=lambda x: -x[1])[0]
    return winner[0]


def chim_score(
        diamond_file_path,
        genes_called=0,
        tax_level='',
        sensitive=False,
        min_mapped_genes=11,
        db="progenomes_2.1",
):
    diamond_df = read_diamond_output(diamond_file_path)
    base_data = create_base_data(diamond_df, db)
    genes_mapped, contig_count = get_stats(diamond_df)
    genome_name = os.path.basename(diamond_file_path).split(".diamond.")[0]
    abundant_lineages_cutoff = get_abundant_lineages_cutoff(sensitive, genes_mapped)

    score = get_scores_for_taxlevel(
        base_data,
        tax_level,
        abundant_lineages_cutoff,
        genome_name,
        genes_called,
        genes_mapped,
        contig_count,
        min_mapped_genes,
    )

    return

def contig_to_n_genes(path_faa):
    d_faa = read_fasta(path_faa)
    out = defaultdict(lambda: 0)
    for seq_id in d_faa.keys():
        contig_id = seq_id[0:seq_id.rfind('_')]
        out[contig_id] += 1
    return out


def worker(job):
    gid, gtdb_taxonomy, tax_level, db = job
    gid_root = get_gid_root(gid)

    diamond_file_path = os.path.join(gid_root, 'gunc', 'diamond_output', f'{gid}.diamond.{db}.out')

    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_contig_to_gene_count = contig_to_n_genes(path_faa)
    total_gene_count = sum(d_contig_to_gene_count.values())


    # Initialise the data
    diamond_df = read_diamond_output(diamond_file_path)
    base_data = create_base_data(diamond_df, db)


    # Get the base case
    genes_mapped, contig_count = get_stats(diamond_df)
    abundant_lineages_cutoff = get_abundant_lineages_cutoff(False, genes_mapped)
    base_score = get_scores_for_taxlevel(
        base_data,
        tax_level,
        abundant_lineages_cutoff,
        gid,
        total_gene_count,
        genes_mapped,
        contig_count,
        11,
    )
    base_max_css = base_score['clade_separation_score']

    # Get the contig mapping prop
    expected_taxon = gtdb_taxonomy.split(';')[TAX_LEVELS[base_score['taxonomic_level']]]
    path_contig_ass = os.path.join(gid_root, 'gunc', 'gunc_output', f'{gid}.contig_assignments.tsv')
    d_map_prop = get_contig_mapping_proportion_for_taxon(expected_taxon, path_contig_ass, tax_level)

    rows = list()

    # candidates = [frozenset((x, )) for x in d_contig_to_gene_count.keys()]
    for cur_contig, cur_cnt in tqdm(d_contig_to_gene_count.items()):

        # Find contigs that improve the max css
        diamond_df_job = diamond_df[diamond_df['contig'] != cur_contig]
        base_data_job = base_data[base_data['contig'] != cur_contig]

        genes_mapped, contig_count = get_stats(diamond_df_job)
        abundant_lineages_cutoff = get_abundant_lineages_cutoff(False, genes_mapped)

        score = get_scores_for_taxlevel(
            base_data_job,
            tax_level,
            abundant_lineages_cutoff,
            gid,
            total_gene_count - cur_cnt,
            genes_mapped,
            contig_count,
            11,
            )

        score['removed'] = cur_contig
        score['prop'] = d_map_prop[cur_contig]
        score['delta'] = base_max_css - score['clade_separation_score']
        rows.append(dict(score))

    # new_candidates = list()
    # for comb in combinations(decreased_css, 2):
    #     new_candidates.append(comb[0] | comb[1])
    # candidates = new_candidates

    # Generate the queue
    df = pd.DataFrame(rows)
    return

    queue = set()
    keys = tuple(sorted(set_contigs))
    for i in range(1, min(len(keys), 3)):
        for comb in combinations(keys, i):
            queue.add(frozenset(comb))
    out = list()

    # Subset the dataframe to contain only the contigs in the queue
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap_unordered(worker_2, queue), total=len(queue)))


    for job in tqdm(queue):
        job = frozenset(job)

        diamond_df_job = diamond_df[~diamond_df['contig'].isin(job)]
        base_data_job = base_data[~base_data['contig'].isin(job)]

        genes_mapped, contig_count = get_stats(diamond_df_job)
        abundant_lineages_cutoff = get_abundant_lineages_cutoff(False, genes_mapped)

        score = get_scores_for_taxlevel(
            base_data_job,
            tax_level,
            abundant_lineages_cutoff,
            gid,
            total_gene_count - sum([d_contig_to_gene_count[x] for x in job]),
            genes_mapped,
            contig_count,
            11,
        )
        score['removed'] = ','.join(job)
        out.append(dict(score))

    df = pd.DataFrame(out)


    return


def worker_2(task):
    job, diamond_df, base_data, tax_level, gid, total_gene_count, d_contig_to_gene_count = task

    out = list()
    diamond_df_job = diamond_df[~diamond_df['contig'].isin(job)]
    base_data_job = base_data[~base_data['contig'].isin(job)]

    genes_mapped, contig_count = get_stats(diamond_df_job)
    abundant_lineages_cutoff = get_abundant_lineages_cutoff(False, genes_mapped)

    score = get_scores_for_taxlevel(
        base_data_job,
        tax_level,
        abundant_lineages_cutoff,
        gid,
        total_gene_count - sum([d_contig_to_gene_count[x] for x in job]),
        genes_mapped,
        contig_count,
        11,
    )
    score['removed'] = ','.join(job)
    out.append(dict(score))
    return out




def main():
    df = get_df()

    queue = list()
    for gid, row in df.iterrows():
        queue.append((gid, row['gtdb_taxonomy'], row['taxonomic_level'], 'gtdb_95'))

    # with mp.Pool(processes=mp.cpu_count()) as pool:
    #     results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue), smoothing=0.05))

    for job in queue:
        # if job[0] != 'GCA_002729355.1':
        #     continue

        worker(job)
    return


if __name__ == '__main__':
    main()
