import multiprocessing as mp
import os

import pandas as pd
from Bio import SeqIO
from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_MARKER, R207_MARKERS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class GetMarkersOnContigsForFailed(LuigiTask):
    """Get the muq/mul/mis/unq markers for each genome at cutoff values."""

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'markers_on_contigs_for_failed.h5'))

    def run(self):
        log('Getting markers on contigs for GUNC failed genomes', title=True)
        self.make_output_dirs()
        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        queue = list()
        for gid, row in df_merged.iterrows():
            queue.append([
                gid,
                row['source'],
                row['domain'],
                row['taxonomic_level']
            ])

            if DEBUG and len(queue) > 3:
                break

        if DEBUG:
            results = [run_on_gid(x) for x in queue]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(run_on_gid, queue), total=len(queue)))

        log('Creating dataframe')
        rows = list()
        for result in results:
            rows.extend(result)
        df = pd.DataFrame(rows)

        log('Saving dataframe')
        if not DEBUG:
            self.save_hdf(df)


def run_on_gid(job):
    gid, gunc_source, domain, max_css = job

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
    d_pct_to_contigs_to_remove = contigs_to_remove_from_gunc(d_fna, df_contig_assign, taxon, tax_level)

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_faa = dict(read_fasta(path_faa))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # Iterate over each percentage to determine which makers are kept
    base_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th)
    out = [{
        'gid': gid,
        'pct': 0
    }]
    for marker in R207_MARKERS:
        if marker in base_markers['unq']:
            out[0][marker] = 'unq'
        elif marker in base_markers['mul']:
            out[0][marker] = 'mul'
        elif marker in base_markers['muq']:
            out[0][marker] = 'muq'
        elif marker in base_markers['mis']:
            out[0][marker] = 'mis'
        else:
            raise Exception('Marker not found')

    for cur_pct, contigs_to_omit in d_pct_to_contigs_to_remove.items():
        cur_results_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th, contigs_to_omit)
        cur_result = {
            'gid': gid,
            'pct': cur_pct
        }
        for marker in R207_MARKERS:
            if marker in cur_results_markers['unq']:
                cur_result[marker] = 'unq'
            elif marker in cur_results_markers['mul']:
                cur_result[marker] = 'mul'
            elif marker in cur_results_markers['muq']:
                cur_result[marker] = 'muq'
            elif marker in cur_results_markers['mis']:
                cur_result[marker] = 'mis'
            else:
                raise Exception('Marker not found')
        out.append(cur_result)

    return out
