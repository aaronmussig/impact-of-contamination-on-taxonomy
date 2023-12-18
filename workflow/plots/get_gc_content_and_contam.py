import multiprocessing as mp
import os

import pandas as pd
from Bio import SeqIO
from luigi import LocalTarget
from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_PLOTS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, \
    contigs_to_remove_from_gunc_only_pct
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def get_gc_worker(job):
    gid, source, domain, tax_level = job
    return get_contig_contam_prop(gid, source, domain, tax_level)


class GetGcContentAndContam(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_PLOTS, 'gc_and_contig_contam.h5'))

    def run(self):
        log('Getting GC content and contamination% for each contig', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        gid_subset = {'GCA_002170165.1', 'GCA_002401495.1', 'GCA_014385135.2', 'GCA_017515185.1', 'GCA_900765805.1',
                      'GCA_905214645.1'}

        log('Creating queue')
        queue = list()
        for gid, row in tqdm(df_merged.iterrows(), total=len(df_merged)):

            if DEBUG and gid not in gid_subset:
                continue

            queue.append((gid, row['source'], row['domain'], row['taxonomic_level']))

        if DEBUG:
            all_rows = list()
            [all_rows.extend(get_gc_worker(x)) for x in queue]
            df = pd.DataFrame(all_rows)
            return


        else:
            log('Processing queue')
            all_rows = list()
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(get_gc_worker, queue), total=len(queue)))
                for result in results:
                    all_rows.extend(result)

            log('Creating dataframe')
            df = pd.DataFrame(all_rows)

            log('Saving...')
            if not DEBUG:
                self.save_hdf(df)

        return


def get_contig_contam_prop(gid, gunc_source, domain, max_css):
    # Read the FASTA file
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')

    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    if gunc_source == 'gtdb':
        path_contig_assign = os.path.join(gid_root, 'gunc_r95/gunc_output', f'{gid}.contig_assignments.tsv')
        gunc_domain = domain
    elif gunc_source == 'progenomes':
        path_contig_assign = os.path.join(gid_root, 'gunc_pro/gunc_output', f'{gid}.contig_assignments.tsv')
        if domain == 'd__Bacteria':
            gunc_domain = '2 Bacteria'
        elif domain == 'd__Archaea':
            gunc_domain = '2157 Archaea'
        else:
            raise ValueError(f'Unknown domain: {domain}')
    else:
        raise Exception(f'Unknown gunc source: {gunc_source}')
    df_contig_assign = read_contig_assignments_tsv(path_contig_assign)

    # Determine the percentage values at which this genome can have contigs removed
    taxon, tax_level = get_taxonomy_by_majority_vote_gunc(df_contig_assign, max_css, gunc_domain)
    rows = contigs_to_remove_from_gunc_only_pct(gid, d_fna, df_contig_assign, taxon, tax_level)
    return rows
