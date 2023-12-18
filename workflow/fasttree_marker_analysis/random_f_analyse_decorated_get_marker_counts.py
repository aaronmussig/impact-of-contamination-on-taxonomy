import multiprocessing as mp
import os
from collections import defaultdict

import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, FASTTREE_MARKER_ANALYSIS_N_BATCHES
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_marker_analysis.random_d_analyse_decorated import FastTreeMarkerAnalyseDecoratedRandom
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastTreeMarkerAnalyseDecoratedRandomGetMarkerCounts(LuigiTask):
    target_pct = luigi.FloatParameter()

    def requires(self):
        out = {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

        for i in range(FASTTREE_MARKER_ANALYSIS_N_BATCHES):
            out[f'analysis_{i}'] = FastTreeMarkerAnalyseDecoratedRandom(batch_id=i, target_pct=self.target_pct)
        return out

    def output(self):
        return LocalTargetHdf5(
            os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, 'random_f_analyse_decorated_get_marker_counts.h5'))

    def run(self):
        log(f'Analysing results of decorated tree (random) (marker counts)  (pct={self.target_pct})',
            title=True)
        self.make_output_dirs()

        log('Loading previous results')
        d_gid_to_batch_and_set = defaultdict(dict)
        for batch_id in range(1 if DEBUG else FASTTREE_MARKER_ANALYSIS_N_BATCHES):

            log(f'Processing batch id={batch_id}')
            df = pd.read_csv(self.input()[f'analysis_{batch_id}'].path, sep='\t')
            df.set_index('gid', inplace=True)

            # log('Breaking into sets')
            gids_congruent = set()
            gids_incongruent = set()
            gids_no_markers = set()

            # Congruent
            gids_congruent.update(set(df[df['classification'] == 'same_msa'].index))
            gids_congruent.update(set(df[(df['classification'] == 'run_on') & (df['tax_result'] == 'correct')].index))
            gids_congruent.update(set(df[(df['classification'] == 'run_on') & (df['tax_result'] == 'congruent')].index))

            # Incongruent
            gids_incongruent.update(set(df[(df['classification'] == 'run_on') & (df['tax_result'] == 'check')].index))
            gids_incongruent.update(set(df[df['classification'] == 'changed_domain'].index))

            # No markers
            gids_no_markers.update(set(df[df['classification'] == 'no_markers'].index))

            # log('Sanity check')
            assert (gids_congruent.union(gids_no_markers).union(gids_incongruent) == set(df.index))
            assert (len(gids_congruent.intersection(gids_no_markers).intersection(gids_incongruent)) == 0)

            for gid in gids_congruent:
                d_gid_to_batch_and_set[gid][batch_id] = 'congruent'
            for gid in gids_incongruent:
                d_gid_to_batch_and_set[gid][batch_id] = 'incongruent'
            for gid in gids_no_markers:
                d_gid_to_batch_and_set[gid][batch_id] = 'no_markers'

        log('Collecting results')
        queue = sorted(list(d_gid_to_batch_and_set.keys()))
        if DEBUG:
            queue = queue[:55]
        d_gid_to_base_markers = dict()
        with mp.Pool(processes=mp.cpu_count()) as pool:
            for result in tqdm(pool.imap_unordered(get_base_markers_for_gid, queue), total=len(queue)):
                d_gid_to_base_markers[result[0]] = result[1]

        log('Calculating data')
        out = list()
        for batch_id in range(1 if DEBUG else FASTTREE_MARKER_ANALYSIS_N_BATCHES):
            log(f'Processing batch id={batch_id}')
            df = pd.read_csv(self.input()[f'analysis_{batch_id}'].path, sep='\t')
            df.set_index('gid', inplace=True)

            for gid, row in tqdm(df.iterrows(), total=len(df)):

                if DEBUG and gid not in d_gid_to_base_markers:
                    break

                base_markers = d_gid_to_base_markers[gid]

                set_lost = set()
                if len(row['markers_lost']) > 5:
                    set_lost = marker_str_to_set(row['markers_lost'])
                set_gained = set()
                if len(row['markers_gained']) > 5:
                    set_gained = marker_str_to_set(row['markers_gained'])

                new_markers = (base_markers - set_lost).union(set_gained)
                jacc = len(base_markers.intersection(new_markers)) / len(base_markers.union(new_markers))

                out.append({
                    'batch_id': batch_id,
                    'set': d_gid_to_batch_and_set[gid][batch_id],
                    'gid': gid,
                    'original_markers': '|'.join(sorted(base_markers)),
                    'new_markers': '|'.join(sorted(new_markers)) if len(new_markers) > 0 else 'N/A',
                    'jaccard': jacc,
                })

        log('Creating dataframe')
        df_out = pd.DataFrame(out)
        df_out.sort_values(by=['batch_id', 'set', 'gid', 'jaccard'], inplace=True, ignore_index=True)

        if not DEBUG:
            self.save_hdf(df_out)

        return


def get_base_markers_for_gid(job):
    gid = job
    genome = Genome(gid)
    d_marker_hits = genome.get_marker_hits()
    out = set()
    out.update(set(d_marker_hits['unq'].keys()))
    out.update(set(d_marker_hits['muq'].keys()))
    return gid, frozenset(out)


def marker_str_to_set(in_str):
    return frozenset(in_str.lstrip("('").rstrip("',)").split(';'))
