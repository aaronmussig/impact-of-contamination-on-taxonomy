import os

import luigi
import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTANI_RANDOM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani_random.c_run_fastani import FastAniRandomRunFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastAniRandomReportResults(LuigiTask):
    target_pct = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        out = {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'sp_clusters': GtdbSpClustersR207(),
        }
        for i in range(10):
            out[f'ani_{i}'] = FastAniRandomRunFastAni(batch_id=i, target_pct=50)
        return out

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                DIR_OUT_FASTANI_RANDOM,
                f'results_report_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}.h5'
            ))

    def load_ani_results_merged(self):
        out = list()
        for i in tqdm(range(10)):
            cur_df = self.input()[f'ani_{i}'].maybe_read_cached()
            cur_df = cur_df[cur_df['ani'] > 0.0]
            cur_df = cur_df[['uid', 'query', 'ref', 'ani', 'af']]
            cur_df['batch_id'] = i
            out.append(cur_df)
        return pd.concat(out, ignore_index=True)

    def run(self):
        log(f'FastAniRandomReportResults(p={self.target_pct},mash_k={self.mash_k}, mash_s={self.mash_s})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading species cluster information')
        df_sp_cluster = self.input()['sp_clusters'].maybe_read_cached()

        log('Loading ANI results')
        all_results = list()
        for batch_id in range(10):

            log(f'Processing batch id: {batch_id}')
            cur_df = self.input()[f'ani_{batch_id}'].maybe_read_cached()
            cur_df = cur_df[cur_df['ani'] > 0.0]
            cur_df = cur_df[['uid', 'query', 'ref', 'ani', 'af']]
            unique_uids = sorted(set(cur_df['uid'].unique()))
            cur_df.rename(columns={'ref': 'reference'}, inplace=True)
            cur_df.set_index('uid', inplace=True)

            log('Getting information for each genome in the batch file')
            results = list()
            for uid in tqdm(unique_uids):
                df_subset = cur_df.loc[[uid]]
                query_gid = df_subset.iloc[0]['query']

                expected_species = df_meta.loc[query_gid, 'species']
                sp_cluster_info = df_sp_cluster.loc[expected_species]

                gid_rep = sp_cluster_info['rep_genome'][3:]
                ani_radius = float(sp_cluster_info['ani_radius'])

                results.append(run_sp_clustering_on_fastani_values(query_gid, gid_rep, df_subset, ani_radius))

            cur_results_df = pd.DataFrame(results)
            cur_results_df['batch_id'] = batch_id
            all_results.append(cur_results_df)

        log('Creating single dataframe')
        df_all = pd.concat(all_results, ignore_index=True)

        if not DEBUG:
            log('Saving')
            self.save_hdf(df_all)
        return
