"""
This file will compress all data in preparation for upload to Zenodo.

A companion script is available that will decompress the data into the expected format.
"""

import os
import tempfile
from typing import Tuple, Collection
import glob

from magna.util.disk import get_file_size_fmt, move_file
from tqdm import tqdm

from workflow.config import DIR_OUT_ROOT, DIR_OUT_BOOTSTRAP, DIR_OUT_EXTERNAL, DIR_OUT_FASTANI, DIR_OUT_FASTANI_INTER, \
    DIR_OUT_GUNC, DIR_OUT_MARKER, DIR_OUT_PPLACER
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.util.log import log
from workflow.util.paths import get_gid_root
import subprocess

DIR_OUT = '/tmp/gunc-export'
os.makedirs(DIR_OUT, exist_ok=True)


# ---------------------------------------------------------------------------- #
# GENOME DIRECTORY
# ---------------------------------------------------------------------------- #

"""
output/genome/gid.fna
output/genome/prodigal/gid.fna
"""
PATH_NUC_ARCHIVE = os.path.join(DIR_OUT, 'genome_nucleotide.tar.gz')

"""
output/genome/prodigal/gid.faa
"""
PATH_AA_ARCHIVE = os.path.join(DIR_OUT, 'genome_aminos.tar.gz')

"""
output/genome/prodigal/gid_pfam_tophit.tsv
output/genome/prodigal/gid_tigrfam_tophit.tsv
"""
PATH_TOPHIT_ARCHIVE = os.path.join(DIR_OUT, 'genome_tophit.tar.gz')

"""
output/genome/gunc_pro/gene_counts.json
output/genome/gunc_pro/GUNC.progenomes_2.1.maxCSS_level.tsv
output/genome/gunc_pro/diamond_output/{gid}.diamond.progenomes_2.1.out
output/genome/gunc_pro/gunc_output/{gid}.contig_assignments.tsv
output/genome/gunc_pro/gunc_output/{gid}.progenomes_2.1.all_levels.tsv

output/genome/gunc_r95/gene_counts.json
output/genome/gunc_r95/GUNC.gtdb_95.maxCSS_level.tsv
output/genome/gunc_r95/diamond_output/{gid}.diamond.gtdb_95.out
output/genome/gunc_r95/gunc_output/{gid}.contig_assignments.tsv
output/genome/gunc_r95/gunc_output/{gid}.gtdb_95.all_levels.tsv
"""
PATH_GUNC_ARCHIVE = os.path.join(DIR_OUT, 'genome_gunc.tar.gz')

"""
output/genome/closest_100_representatives.h5
"""
PATH_GENOME_H5_FILES = os.path.join(DIR_OUT, 'genome_h5_files.tar.gz')

"""
output/genome/*.h5
output/genome/msa_at_cutoff_values.tsv
"""
PATH_GENOME_TOP_LEVEL_FILES = os.path.join(DIR_OUT, 'genome_top_level_files.tar.gz')


# ---------------------------------------------------------------------------- #
# MISC DIRECTORY
# ---------------------------------------------------------------------------- #
PATH_BOOTSTRAP_ARCHIVE = os.path.join(DIR_OUT, 'bootstrap_dir.tar.gz')
PATH_EXTERNAL_ARCHIVE = os.path.join(DIR_OUT, 'external_dir.tar.gz')
PATH_FASTANI_ARCHIVE = os.path.join(DIR_OUT, 'fastani_dir.tar.gz')
PATH_FASTANI_INTER_ARCHIVE = os.path.join(DIR_OUT, 'fastani_interspecies_dir.tar.gz')
PATH_GUNC_DIR_ARCHIVE = os.path.join(DIR_OUT, 'gunc_dir.tar.gz')
PATH_MARKER_ARCHIVE = os.path.join(DIR_OUT, 'marker_dir.tar.gz')
PATH_PPLACER_ARCHIVE = os.path.join(DIR_OUT, 'pplacer_dir.tar.gz')

"""
Utility methods
"""

def compress_file_paths_into_archive(file_list: Collection[str], path_out: str):
    log(f'Compressing {len(file_list):,} files into {path_out}')

    with tempfile.TemporaryDirectory() as tmp_dir:
        path_input = os.path.join(tmp_dir, 'input.txt')
        with open(path_input, 'w') as f:
            for path in file_list:
                f.write(f'{path}\n')

        compress_pigz(path_input, path_out, n_total=len(file_list))
    return

def compress_pigz(file_list, path_out, n_total: int):
    # cmd = f"tar -cvf - --dereference -T {file_list} | pigz --best > {path_out}"
    cmd = f"tar -cvf - --no-acls --owner=root --group=root -T {file_list} --dereference --transform 's,srv/home/uqamussi/projects/gunc-chimeras/output/,output/,' | pigz --best --iterations 100 --maxsplits 100 > {path_out}"

    with tqdm(total=n_total, unit=' file', desc='Compressing', smoothing=0.1) as p_bar:
        with subprocess.Popen(cmd, stderr=subprocess.PIPE, encoding='utf-8', shell=True) as proc:
            while True:
                line = proc.stderr.readline()
                if not line:
                    break
                p_bar.update()
            proc.wait()

    log(f'Done: {get_file_size_fmt(path_out)}')
    return


def save_hdf5(df, path):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, 'out.h5')
        df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')
        log(f'Copying {get_file_size_fmt(tmp_path)} to {path}')
        move_file(tmp_path, path, checksum=True)


# ---------------------------------------------------------------------------- #




"""
Methods for files that reside in the genome/. directory
"""

def compress_nucleotide_files(gids, path_out):
    if os.path.isfile(path_out):
        log(f'File exists, skipping: {path_out}')
        return

    queue = list()
    for gid in gids:
        gid_root = get_gid_root(gid)
        fna_path = os.path.join(gid_root, f'{gid}.fna')
        gene_fna_path = os.path.join(gid_root, 'prodigal', f'{gid}.fna')
        queue.append(fna_path)
        queue.append(gene_fna_path)

    compress_file_paths_into_archive(queue, path_out)
    return

def compress_amino_acid_files(gids, path_out):
    if os.path.isfile(path_out):
        log(f'File exists, skipping: {path_out}')
        return

    queue = list()
    for gid in gids:
        gid_root = get_gid_root(gid)
        gene_faa_path = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
        queue.append(gene_faa_path)

    compress_file_paths_into_archive(queue, path_out)
    return

def compress_tophit_files(gids, path_out):
    if os.path.isfile(path_out):
        log(f'File exists, skipping: {path_out}')
        return

    queue = list()
    for gid in gids:
        gid_root = get_gid_root(gid)
        pfam_path = os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv')
        tigr_path = os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv')

        queue.append(pfam_path)
        queue.append(tigr_path)

    compress_file_paths_into_archive(queue, path_out)
    return

def compress_gunc_archive(gids, path_out):
    if os.path.isfile(path_out):
        log(f'File exists, skipping: {path_out}')
        return

    queue = list()
    for gid in gids:
        gid_root = get_gid_root(gid)

        queue.append(os.path.join(gid_root, 'gunc_pro', 'gene_counts.json'))
        queue.append(os.path.join(gid_root, 'gunc_pro', 'GUNC.progenomes_2.1.maxCSS_level.tsv'))
        queue.append(os.path.join(gid_root, 'gunc_pro', 'diamond_output', f'{gid}.diamond.progenomes_2.1.out'))
        queue.append(os.path.join(gid_root, 'gunc_pro', 'gunc_output', f'{gid}.contig_assignments.tsv'))
        queue.append(os.path.join(gid_root, 'gunc_pro', 'gunc_output', f'{gid}.progenomes_2.1.all_levels.tsv'))

        queue.append(os.path.join(gid_root, 'gunc_r95', 'gene_counts.json'))
        queue.append(os.path.join(gid_root, 'gunc_r95', 'GUNC.gtdb_95.maxCSS_level.tsv'))
        queue.append(os.path.join(gid_root, 'gunc_r95', 'diamond_output', f'{gid}.diamond.gtdb_95.out'))
        queue.append(os.path.join(gid_root, 'gunc_r95', 'gunc_output', f'{gid}.contig_assignments.tsv'))
        queue.append(os.path.join(gid_root, 'gunc_r95', 'gunc_output', f'{gid}.gtdb_95.all_levels.tsv'))

    compress_file_paths_into_archive(queue, path_out)
    return

def compress_genome_h5_files(path_out):
    queue = [
        os.path.join(DIR_OUT_ROOT, 'genome', 'closest_100_representatives.h5')
    ]
    compress_file_paths_into_archive(queue, path_out)
    return



def compress_entire_directory_recursive(dir_in, path_out):
    queue = glob.glob(dir_in + '/**/*', recursive=True)
    compress_file_paths_into_archive(queue, path_out)
    return

def compress_top_level_files_in_failed_genome_dir(gids, path_out):
    queue = list()
    for gid in gids:
        gid_root = get_gid_root(gid)
        queue.extend(glob.glob(gid_root + '/*.h5'))
        queue.append(os.path.join(gid_root, 'msa_at_cutoff_values.tsv'))
    compress_file_paths_into_archive(queue, path_out)
    return



"""
Main method
"""


def main():
    # Load all genomes
    log('Loading genome list')
    gids_all = tuple(sorted(set(GtdbMetadataR207().output().read().index)))
    log(f'Loaded {len(gids_all):,} genomes')

    log('Loading genomes that failed GUNC')
    gids_fail = tuple(sorted(set(AggregateMaxCssLevelMerged().output().read().index)))
    log(f'Loaded {len(gids_fail):,} genomes')

    gids_all = gids_all[0:10]
    gids_fail = gids_fail[0:10]

    # TODO: Running on page
    # log('Compressing nucleotide files...')
    # compress_nucleotide_files(gids_all, PATH_NUC_ARCHIVE)

    # TODO: Running on page
    # log('Compressing amino acid files...')
    # compress_amino_acid_files(gids_all, PATH_AA_ARCHIVE)

    # TODO: Running on page
    # log('Compressing tophit files...')
    # compress_tophit_files(gids_fail, PATH_TOPHIT_ARCHIVE)

    # TODO: Untested below

    log('Compressing top level files in failed genome directories...')
    compress_top_level_files_in_failed_genome_dir(gids_fail, PATH_GENOME_TOP_LEVEL_FILES)

    log('Compressing GUNC archives...')
    compress_gunc_archive(gids_all, PATH_GUNC_ARCHIVE)

    log('Compressing genome.h5 files...')
    compress_genome_h5_files(PATH_GENOME_H5_FILES)

    log('Compressing bootstrap directory...')
    compress_entire_directory_recursive(DIR_OUT_BOOTSTRAP, PATH_BOOTSTRAP_ARCHIVE)

    log('Compressing external directory...')
    compress_entire_directory_recursive(DIR_OUT_EXTERNAL, PATH_EXTERNAL_ARCHIVE)

    log('Compressing FastANI directory...')
    compress_entire_directory_recursive(DIR_OUT_FASTANI, PATH_FASTANI_ARCHIVE)

    log('Compressing FastANI interspecies directory...')
    compress_entire_directory_recursive(DIR_OUT_FASTANI_INTER, PATH_FASTANI_INTER_ARCHIVE)

    log('Compressing gunc directory...')
    compress_entire_directory_recursive(DIR_OUT_GUNC, PATH_GUNC_DIR_ARCHIVE)

    log('Compressing marker directory...')
    compress_entire_directory_recursive(DIR_OUT_MARKER, PATH_MARKER_ARCHIVE)

    log('Compressing pplacer directory...')
    compress_entire_directory_recursive(DIR_OUT_PPLACER, PATH_PPLACER_ARCHIVE)
    return


if __name__ == '__main__':
    main()
