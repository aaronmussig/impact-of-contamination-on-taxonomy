import multiprocessing as mp
import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_MARKER, PCT_VALUES
from workflow.config import R207_AR53_HMM, R207_BAC120_HMM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid
from workflow.method.randomly_select_contigs_up_to_pct import randomly_select_contigs_up_to_pct
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.tophit import TopHitTigrFile, TopHitPfamFile
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root

N_REPEATS = 100


class RandomContigRemovalToMarkersPresent(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'random_contig_removal_to_markers_present.h5'))

    def run(self):
        log('Randomly removing contigs for each genome to determine what markers remain', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        queue = list()
        for gid, row in df_merged.iterrows():
            queue.append([gid, row['domain']])

            if DEBUG and len(queue) > 3:
                break
        log(f'Enqueued {len(queue):,} genomes')

        if DEBUG:
            results = [run_on_gid(x) for x in queue]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(run_on_gid, queue), total=len(queue)))

        log('Converting results into rows')
        rows = list()
        [rows.extend(x) for x in results]

        log('Saving results')
        df = pd.DataFrame(rows)
        df.sort_values(by=['gid', 'pct'], inplace=True)

        if not DEBUG:
            self.save_hdf(df)
        return


def run_on_gid(job):
    gid, domain = job

    # Set the markers for the current domain
    if domain == 'd__Bacteria':
        marker_set = frozenset(R207_BAC120_HMM)
    elif domain == 'd__Archaea':
        marker_set = frozenset(R207_AR53_HMM)
    else:
        raise ValueError(f'Invalid domain: {domain}')

    # Read the FASTA file
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')

    # Load the FNA
    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_faa = dict(read_fasta(path_faa))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # Iterate over each percentage to determine which makers are kept
    base_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th)
    base_markers_set = set({**base_markers['unq'], **base_markers['muq']}.keys())
    base_markers_for_domain = base_markers_set.intersection(marker_set)
    base_markers_pct = 100 * len(base_markers_for_domain) / len(marker_set)
    out = [{
        'gid': gid,
        'pct': 0,
        'pct_markers': base_markers_pct,
    }]

    # Randomly select markers up to X% of the genome removed
    for pct in PCT_VALUES:

        # Repeat
        pct_vals = list()
        for n in range(N_REPEATS):
            # Get the marker count
            contigs_to_keep, contigs_to_remove = randomly_select_contigs_up_to_pct(d_fna, pct)
            cur_results_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th, contigs_to_remove)
            valid_markers = {**cur_results_markers['muq'], **cur_results_markers['unq']}
            valid_domain_markers = marker_set.intersection(frozenset(valid_markers.keys()))

            # Store the percentage
            pct_vals.append(100 * len(valid_domain_markers) / len(marker_set))

        # Store the average
        pct_avg = float(np.mean(pct_vals))
        out.append({
            'gid': gid,
            'pct': pct,
            'pct_markers': pct_avg,
        })

    return out
