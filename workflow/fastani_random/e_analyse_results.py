import os
import tempfile

import luigi
import pandas as pd
from luigi import LocalTarget
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE, DEBUG, DIR_OUT_FASTANI_RANDOM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani_congruence.b_run_fastani import FastAniCongruenceRunFastAni
from workflow.fastani_congruence.c_analyse_fastani import FastAniCongruenceAnalyseFastAni
from workflow.fastani_random.d_report_results import FastAniRandomReportResults
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.genome_pct_congruence_contigs_removed import GenomePctCongruenceContigsRemoved
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from collections import Counter, defaultdict
import numpy as np

from workflow.v2_fastani_rep_to_closest_rep.a_create_jobs import V2FastAniRepToClosestRepCreateJobs


class FastAniRandomAnalyseResults(LuigiTask):
    target_pct = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'sp_results': FastAniRandomReportResults(target_pct=self.target_pct, mash_k=self.mash_k, mash_s=self.mash_s),
            'fastani_rep':  V2FastAniRepToClosestRepCreateJobs()
        }

    def output(self):
        return LocalTarget(
            os.path.join(
                DIR_OUT_FASTANI_RANDOM,
                f'results_analysis_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}.tsv'
            ))

    def run(self):
        log(f'FastAniRandomAnalyseResults(p={self.target_pct}, mash_k={self.mash_k}, mash_s={self.mash_s})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_rep = df_meta['gtdb_genome_representative'].to_dict()

        log('Loading distances from reps to the closest rep')
        d_fastani_rep = get_closest_rep(self.input()['fastani_rep'].maybe_read_cached())
        bac_non_rep_gids = set(df_meta[(df_meta['domain'] == 'd__Bacteria') & (df_meta['gtdb_representative'] == 'f')].index)

        log('Loading species clustering results')
        df_sp_results = self.input()['sp_results'].maybe_read_cached()
        print(f'Original length: {len(df_sp_results):,}')
        df_sp_results = df_sp_results[df_sp_results['gid'].isin(bac_non_rep_gids)]
        print(f'After removing reps: {len(df_sp_results):,}')

        all_no_change = list()
        all_new_sp_cluster = list()
        all_changed_species = list()

        # Additionally, we want to know what the distance is for each case to it's expected rep.
        d_batch_to_correct = dict()
        d_batch_to_changed_sp = dict()
        d_batch_to_new_sp = dict()

        for batch_id in df_sp_results['batch_id'].unique():
            cur_batch_rows = df_sp_results[df_sp_results['batch_id'] == batch_id]

            rows_correct = cur_batch_rows[cur_batch_rows['same'] == True]
            rows_wrong = cur_batch_rows[cur_batch_rows['same'] == False]
            rows_changed_sp = rows_wrong[rows_wrong['type'] == 'sp_rep']
            rows_new_sp_cluster = rows_wrong[rows_wrong['type'].isin({'no_ani', 'no_af'})]

            gids_correct = set(rows_correct['gid'])
            gids_changed_sp = set(rows_changed_sp['gid'])
            gids_new_sp_cluster = set(rows_new_sp_cluster['gid'])

            cur_n_no_change = len(rows_correct)
            cur_n_new_sp_cluster = len(rows_new_sp_cluster)
            cur_n_changed_species = len(rows_changed_sp)

            all_no_change.append(cur_n_no_change / len(cur_batch_rows) * 100)
            all_new_sp_cluster.append(cur_n_new_sp_cluster / len(cur_batch_rows) * 100)
            all_changed_species.append(cur_n_changed_species / len(cur_batch_rows) * 100)

            # For each of the genomes collect the representative genome ANI
            d_batch_to_correct[int(batch_id)] = get_rep_and_closest_ani(d_fastani_rep, gids_correct, d_gid_to_rep)
            d_batch_to_changed_sp[int(batch_id)] = get_rep_and_closest_ani(d_fastani_rep, gids_changed_sp, d_gid_to_rep)
            d_batch_to_new_sp[int(batch_id)] = get_rep_and_closest_ani(d_fastani_rep, gids_new_sp_cluster, d_gid_to_rep)

        no_change_mean = np.mean(all_no_change)
        no_change_std = np.std(all_no_change)

        new_sp_cluster_mean = np.mean(all_new_sp_cluster)
        new_sp_cluster_std = np.std(all_new_sp_cluster)

        changed_species_mean = np.mean(all_changed_species)
        changed_species_std = np.std(all_changed_species)

        log(f'No change: {no_change_mean:.2f} ± {no_change_std:.2f}')
        log(f'New sp cluster: {new_sp_cluster_mean:.2f} ± {new_sp_cluster_std:.2f}')
        log(f'Changed species: {changed_species_mean:.2f} ± {changed_species_std:.2f}')


        correct_ani_vals = [np.mean(d_batch_to_correct[x]) for x in range(10)]
        changed_sp_ani_vals = [np.mean(d_batch_to_changed_sp[x]) for x in range(10)]
        new_sp_ani_vals = [np.mean(d_batch_to_new_sp[x]) for x in range(10)]

        correct_ani_mean = np.mean(correct_ani_vals)
        correct_ani_std = np.std(correct_ani_vals)

        changed_sp_ani_mean = np.mean(changed_sp_ani_vals)
        changed_sp_ani_std = np.std(changed_sp_ani_vals)

        new_sp_ani_mean = np.mean(new_sp_ani_vals)
        new_sp_ani_std = np.std(new_sp_ani_vals)

        log(f'Correct ANI: {100-correct_ani_mean:.2f} ± {correct_ani_std:.2f}')
        log(f'Changed sp ANI: {100-changed_sp_ani_mean:.2f} ± {changed_sp_ani_std:.2f}')
        log(f'New sp ANI: {100-new_sp_ani_mean:.2f} ± {new_sp_ani_std:.2f}')

        return

def get_closest_rep(df):
    out = dict()
    for row in tqdm(df.itertuples(), total=len(df)):
        if row.af < 0.5 or row.query == row.ref:
            continue

        if row.query in out:
            if row.ani > out[row.query][1]:
                out[row.query] = (row.ref, row.ani)
        else:
            out[row.query] = (row.ref, row.ani, row.af)

        if row.ref in out:
            if row.ani > out[row.ref][1]:
                out[row.ref] = (row.query, row.ani)
        else:
            out[row.ref] = (row.query, row.ani)
    return out


def get_rep_and_closest_ani(d_rep_ani, gids, d_gid_to_rep):
    out = list()
    n_skipped = 0
    for gid in gids:
        rep = d_gid_to_rep[gid][3:]
        hit = d_rep_ani.get(rep)
        if hit is None:
            n_skipped += 1
            continue
        out.append(hit[1])

    print(f'Skipped: {n_skipped:,}')
    return out
