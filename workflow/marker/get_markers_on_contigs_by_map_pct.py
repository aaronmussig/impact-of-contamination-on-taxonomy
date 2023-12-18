import multiprocessing as mp
import os
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_MARKER, R207_MARKERS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc, \
    contigs_to_remove_from_gunc_only_pct
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid, get_marker_hits_for_gid_include_mul
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class GetMarkersOnContigsByMapPct(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'markers_on_contigs_for_failed_by_map_pct_v2.h5'))

    def run(self):
        log('Getting markers on contigs for GUNC failed genomes (by map pct)', title=True)
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

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_faa = dict(read_fasta(path_faa))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # Iterate over each percentage to determine which makers are kept
    base_marker, marker_df = get_marker_hits_for_gid_include_mul(d_faa, pfam_th, tigr_th)
    out = list()
    for cur_row in marker_df:
        cur_gene_id = cur_row['hit'].gene_id
        cur_contig = cur_gene_id[:cur_gene_id.rindex('_')]
        cur_marker = cur_row['hit'].hmm_id
        cur_marker_n = cur_row['n']

        if cur_marker in base_marker['unq']:
            marker_type = 'unq'
        elif cur_marker in base_marker['mul']:
            marker_type = 'mul'
        elif cur_marker in base_marker['muq']:
            marker_type = 'muq'
        else:
            raise ValueError(f'Unknown marker type for {cur_marker}')

        out.append({
            'gid': gid,
            'marker': cur_marker,
            'n': cur_marker_n,
            'marker_type': marker_type,
            'contig': cur_contig,
        })

    return out
