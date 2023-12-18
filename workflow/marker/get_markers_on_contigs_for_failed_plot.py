import multiprocessing as mp
import os
from collections import defaultdict, Counter

import pandas as pd
from Bio import SeqIO
from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_MARKER, R207_MARKERS, PCT_VALUES
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_markers_on_contigs_for_failed import GetMarkersOnContigsForFailed
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root
import seaborn as sns
import matplotlib.pyplot as plt


class GetMarkersOnContigsForFailedPlot(LuigiTask):
    """Get the muq/mul/mis/unq markers for each genome at cutoff values."""

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'marker_info': GetMarkersOnContigsForFailed(),
        }

    def output(self):
        return None

    def run(self):
        log('Getting markers on contigs for GUNC failed genomes (plotting)', title=True)

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        log('Loading results from previous step')
        df_marker = self.input()['marker_info'].read() if not DEBUG else self.input()['marker_info'].read_cached()

        log('Converting to dictionary')
        gid_to_rows = defaultdict(dict)
        for _, row in tqdm(df_marker.iterrows(), total=len(df_marker)):
            gid_to_rows[row['gid']][row['pct']] = tuple([row[x] for x in R207_MARKERS])

        log('Find those rows that differ from the baseline')
        diff_rows = get_rows_that_differ_from_baseline(gid_to_rows)

        plot_total_number_of_markers_at_each_pct_value(gid_to_rows)

        log('Load the baseline data')



        return

def get_rows_that_differ_from_baseline(gid_to_rows):
    differ_rows = list()

    for gid, d_pct_to_row in tqdm(gid_to_rows.items()):
        first_row = d_pct_to_row[0]
        prev_row = first_row

        for pct in PCT_VALUES:
            cur_row = d_pct_to_row.get(pct, None)

            # Not present, so carry forward the previous value
            if cur_row is None:
                use_row = prev_row

            # Is present, so use the current value
            else:
                use_row = cur_row
                prev_row = use_row

            if use_row != first_row:
                differ_rows.append((gid, pct))
    return differ_rows


def plot_total_number_of_markers_at_each_pct_value(gid_to_rows):
    # Note: Carry forward the previous cutoff value data if no further data are present

    d_pct_to_counts = defaultdict(lambda: defaultdict(lambda: 0))

    for gid, d_pct_to_row in tqdm(gid_to_rows.items()):

        prev_count = Counter(d_pct_to_row[0])
        d_pct_to_counts[0]['muq'] += prev_count.get('muq', 0)
        d_pct_to_counts[0]['mis'] += prev_count.get('mis', 0)
        d_pct_to_counts[0]['mul'] += prev_count.get('mul', 0)
        d_pct_to_counts[0]['unq'] += prev_count.get('unq', 0)

        for pct in PCT_VALUES:
            cur_row = d_pct_to_row.get(pct, None)

            # Not present, so carry forward the previous value
            if cur_row is None:
                use_count = prev_count

            # Is present, so use the current value
            else:
                use_count = Counter(cur_row)
                prev_count = use_count

            d_pct_to_counts[pct]['muq'] += use_count.get('muq', 0)
            d_pct_to_counts[pct]['mis'] += use_count.get('mis', 0)
            d_pct_to_counts[pct]['mul'] += use_count.get('mul', 0)
            d_pct_to_counts[pct]['unq'] += use_count.get('unq', 0)

    rows = list()
    for pct, d_cnts in d_pct_to_counts.items():
        for marker_type, marker_cnt in d_cnts.items():
            rows.append({
                'pct': pct,
                'marker_type': marker_type,
                'marker_cnt': marker_cnt,
            })
    df = pd.DataFrame(rows)

    sns.lineplot(x='pct', y='marker_cnt', hue='marker_type', data=df)
    plt.show()

    return

"""
Things we're trying to do here:


For this one:
- gid breakdown at pct, what prop are reps?
- how many in total?



1. At each pct, can we still infer the same domain?

2a. For reps: At each pct, do we still have the genome the same spot in the tree?
2b. For other: At each pct, can we still place it back in the tree


Things to do:
1. Generate the MSA at a pct value
2. For reps, infer the tree without that genome in it (lots of trees!)
3. For non-reps, use gtdbtk to place it into the tree

 

"""
