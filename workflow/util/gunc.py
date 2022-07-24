import json
import os
from collections import defaultdict
from typing import Dict, Literal

import numpy as np
import pandas as pd
from magna.util.disk import copy_file

from workflow.config import TAX_LEVELS
from workflow.util.log import log


def move_gene_counts(path: str, gid_target_dirs: Dict[str, str]):
    with open(path, 'r') as f:
        data = json.load(f)
    for gid, gene_count in data.items():
        target_dir = gid_target_dirs[gid]
        new_path = os.path.join(target_dir, 'gene_counts.json')
        with open(new_path, 'w') as f:
            json.dump({gid: gene_count}, f, indent=4)


def move_diamond_output(diamond_dir, gid_target_dirs: Dict[str, str], db):
    for gid, target_dir in gid_target_dirs.items():
        diamond_path = os.path.join(diamond_dir, f'{gid}.diamond.{db}.out')
        new_path = os.path.join(target_dir, 'diamond_output', f'{gid}.diamond.{db}.out')
        copy_file(diamond_path, new_path, checksum=True)


def move_maxcss_level(path, gid_target_dirs: Dict[str, str], db):
    gids_done = set()
    with open(path, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            cols = line.strip().split('\t')
            gid = cols[0]
            new_path = os.path.join(gid_target_dirs[gid], f'GUNC.{db}.maxCSS_level.tsv')
            with open(new_path, 'w') as fw:
                fw.write(header)
                fw.write('\t'.join(cols) + '\n')
            gids_done.add(gid)
    if gids_done != set(gid_target_dirs.keys()):
        raise Exception(f'Not all gids were found in {path}')


def move_contig_assignments(gunc_dir, gid_target_dirs: Dict[str, str]):
    for gid, target_dir in gid_target_dirs.items():
        path_contig = os.path.join(gunc_dir, f'{gid}.contig_assignments.tsv')
        path_contig_new = os.path.join(target_dir, 'gunc_output', f'{gid}.contig_assignments.tsv')
        copy_file(path_contig, path_contig_new, checksum=True)
    return


def move_all_levels(gunc_dir, gid_target_dirs: Dict[str, str], db):
    for gid, target_dir in gid_target_dirs.items():
        path_all_levels = os.path.join(gunc_dir, f'{gid}.{db}.all_levels.tsv')
        path_all_levels_new = os.path.join(target_dir, 'gunc_output', f'{gid}.{db}.all_levels.tsv')
        copy_file(path_all_levels, path_all_levels_new, checksum=True)
    return


def create_output_dirs(gid_target_dirs: Dict[str, str]):
    log('Creating output directories...')
    target_dirs = set(gid_target_dirs.values())
    for target_dir in target_dirs:
        diamond_dir = os.path.join(target_dir, 'diamond_output')
        gunc_dir = os.path.join(target_dir, 'gunc_output')
        os.makedirs(diamond_dir, exist_ok=True)
        os.makedirs(gunc_dir, exist_ok=True)


def move_gunc_directory(from_dir: str, gid_target_dirs: Dict[str, str], db='gtdb_95'):
    """Parse the result of a GUNC output batch and move to the appropriate location"""
    # Ensure directories exist
    create_output_dirs(gid_target_dirs)

    # Move files
    move_gene_counts(os.path.join(from_dir, 'gene_counts.json'), gid_target_dirs)
    move_diamond_output(os.path.join(from_dir, 'diamond_output'), gid_target_dirs, db=db)
    move_maxcss_level(os.path.join(from_dir, f'GUNC.{db}.maxCSS_level.tsv'), gid_target_dirs, db)
    move_contig_assignments(os.path.join(from_dir, 'gunc_output'), gid_target_dirs)
    move_all_levels(os.path.join(from_dir, 'gunc_output'), gid_target_dirs, db)
    return

def gunc_are_all_files_present(gid: str, gunc_dir: str, db: Literal['gtdb_95', 'progenomes_2.1']):
    """Check if the GUNC directory is valid for a specific gid"""
    gene_counts_path = get_gene_counts_path(gunc_dir)
    diamond_path = get_diamond_path(gunc_dir, gid, db)
    max_css_path = get_max_css_path(gunc_dir, db)
    contig_assign_path = get_contig_assignments_path(gunc_dir, gid)
    all_levels_path = get_all_levels_path(gunc_dir, gid, db)

    for path in [gene_counts_path, diamond_path, max_css_path, contig_assign_path, all_levels_path]:
        try:
            if os.path.getsize(path) < 5:
                return False
        except:
            return False
    return True



def get_max_css_path(gunc_dir: str, db='gtdb_95'):
    return os.path.join(gunc_dir, f'GUNC.{db}.maxCSS_level.tsv')


def get_diamond_path(gunc_dir: str, gid: str, db = 'gtdb_95'):
    return os.path.join(gunc_dir, 'diamond_output', f'{gid}.diamond.{db}.out')


def get_all_levels_path(gunc_dir: str, gid: str, db='gtdb_95'):
    return os.path.join(gunc_dir, 'gunc_output', f'{gid}.{db}.all_levels.tsv')


def get_contig_assignments_path(gunc_dir: str, gid: str):
    return os.path.join(gunc_dir, 'gunc_output', f'{gid}.contig_assignments.tsv')


def get_gene_counts_path(gunc_dir: str):
    return os.path.join(gunc_dir, 'gene_counts.json')

def try_get_taxonomy_from_contig_mapping_file_for_progenomes(contig_assignment_path: str, taxonomic_level: str):
    df_contig = pd.read_csv(contig_assignment_path, sep='\t')
    df_contig = df_contig[df_contig['tax_level'] == taxonomic_level]

    d_tax_to_count = defaultdict(lambda: 0)
    for _, row in df_contig.iterrows():
        d_tax_to_count[row['assignment']] += row['count_of_genes_assigned']

    sorted_counts = sorted(d_tax_to_count.items(), key=lambda x: x[1], reverse=True)

    highest_tax = sorted_counts[0][0]
    highest_cnt = sorted_counts[0][1]

    total_assignments = sum(d_tax_to_count.values())
    cur_map_prop = highest_cnt / total_assignments

    if cur_map_prop <= 0.5:
        raise Exception(f'Ambiguous taxonomy assignment: {contig_assignment_path}')

    return highest_tax

def get_contig_mapping_proportion_for_taxon(expected_taxon: str, contig_assignment_path: str, tax_level: str) -> Dict[str, float]:
    """Calculates the proportion of TRUE taxa that are assigned to a contig at a rank.

    Args:
        expected_taxon: The true taxon of the genome.
        contig_assignment_path: The path to the contig assignment file (h5).
        tax_level: The taxonomic level to check

    Returns:
        The proportion of true mappings for each contig at each taxonomic level.
        e.g. {'contig_a': 0.2, 'contig_b': 0.8, ...}
             0.2% of the true taxa mapped to contig_a
             0.8% of the true taxa mapped to contig_b
    """
    df_contig = pd.read_csv(contig_assignment_path, sep='\t')
    df_contig = df_contig[df_contig['tax_level'] == tax_level]

    d_assignment_counts = defaultdict(list)
    for _, row in df_contig.iterrows():
        d_assignment_counts[row['contig']].append(
            (row['assignment'], row['count_of_genes_assigned'])
        )

    out = dict()
    for contig, assignments in d_assignment_counts.items():
        # Find the proportion of this contig
        total = 0
        n_true = 0
        for taxon, count in assignments:
            if taxon == expected_taxon:
                n_true += count
            total += count

        # Return the proportion of this contig
        out[contig] = n_true / total
    return out


def get_contig_mapping_proportion(taxonomy: str, contig_assignment_path: str) -> Dict[str, Dict[str, float]]:
    """Calculates the proportion of TRUE taxa that are assigned to a contig at each rank.

    Args:
        taxonomy: The true taxonomy of the genome.
        contig_assignment_path: The path to the contig assignment file (h5).

    Returns:
        The proportion of true mappings for each contig at each taxonomic level.
        e.g. {'kingdom': { 'contig_a': 0.2, 'contig_b': 0.8}, ...}
             0.2% of the true taxa mapped to contig_a
             0.8% of the true taxa mapped to contig_b
    """
    df_contig = pd.read_csv(contig_assignment_path, sep='\t')

    d_tax_level_to_assignment_counts = defaultdict(lambda: defaultdict(list))
    for _, row in df_contig.iterrows():
        tax_level = row['tax_level']
        d_tax_level_to_assignment_counts[tax_level][row['contig']].append(
            (row['assignment'], row['count_of_genes_assigned'])
        )

    out = dict()
    expected_taxonomy = taxonomy.split(';')
    for tax_level, d_contigs in d_tax_level_to_assignment_counts.items():
        expected_taxon = expected_taxonomy[TAX_LEVELS[tax_level]]
        out[tax_level] = dict()
        for contig, assignments in d_contigs.items():

            # Find the proportion of this contig
            total = 0
            n_true = 0
            for taxon, count in assignments:
                if taxon == expected_taxon:
                    n_true += count
                total += count

            # Return the proportion of this contig
            out[tax_level][contig] = n_true / total
    return out


MAX_CSS_LEVEL_DTYPE = {
    'genome': object,
    'n_genes_called': np.uintc,
    'n_genes_mapped': np.uintc,
    'n_contigs': np.uintc,
    'taxonomic_level': object,
    'proportion_genes_retained_in_major_clades': np.float64,
    'genes_retained_index': np.float64,
    'clade_separation_score': np.float64,
    'contamination_portion': np.float64,
    'n_effective_surplus_clades': np.float64,
    'mean_hit_identity': np.float64,
    'reference_representation_score': np.float64,
    'pass.GUNC': bool
}

# def main():
#
#     a = '/tmp/pro'
#     b = {'GCA_000006155.2': '/tmp/asd', 'GCA_000007185.1': '/tmp/dsa'}
#     c = 'progenomes_2.1'
#     move_gunc_directory(a, b, c)
#
# if __name__ == '__main__':
#     main()
