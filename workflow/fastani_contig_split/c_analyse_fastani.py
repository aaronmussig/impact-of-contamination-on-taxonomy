import os

import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONTIG_SPLIT
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani_contig_split.b_run_fastani import FastAniContigSplitRunFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastAniContigSplitAnalyseFastAni(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    @property
    def root_dir(self):
        return DIR_OUT_FASTANI_CONTIG_SPLIT

    def requires(self):
        return {
            'ani': FastAniContigSplitRunFastAni(target_pct=self.target_pct, mash_k=self.mash_k, mash_s=self.mash_s),
            'meta': GtdbMetadataR207(),
            'sp_clusters': GtdbSpClustersR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                self.root_dir, f'results_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}.h5'
            ))

    def run(self):
        log(f'FastAniContigSplitAnalyseFastAni(p={self.target_pct}, mash_k={self.mash_k}, mash_s={self.mash_s})',
            title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading species cluster information')
        df_sp_cluster = self.input()['sp_clusters'].maybe_read_cached()

        log('Loading ANI scores')
        df_ani = self.input()['ani'].maybe_read_cached()
        df_ani.rename(columns={'ref': 'reference'}, inplace=True)
        df_ani.set_index('query', inplace=True)

        log('Loading MaxCSS scores')
        df_css = self.input()['max_css'].maybe_read_cached()

        log('Processing failed genomes')
        results = list()
        for gid in tqdm(df_css.index, total=len(df_css.index)):

            for suffix in ('', '_C'):
                gid_suffix = f'{gid}{suffix}'

                expected_species = df_meta.loc[gid, 'species']
                sp_cluster_info = df_sp_cluster.loc[expected_species]

                gid_rep = sp_cluster_info['rep_genome'][3:]

                try:
                    df_fastani_subset = df_ani.loc[[gid_suffix]]
                except KeyError:
                    df_fastani_subset = pd.DataFrame([{
                        'query': gid_suffix,
                        'reference': 'NA',
                        'ani': 0,
                        'af': 0,
                        'ani_qvr': 0,
                        'ani_rvq': 0,
                        'af_qvr': 0,
                        'af_rvq': 0,
                    }])
                    # # This happens if there were no hits within any radius
                    # results.append({
                    #     'gid': gid,
                    #     'new_sp_rep': None,
                    #     'ani': None,
                    #     'af': None,
                    #     'type': 'no_hits',
                    #     'same': True
                    # })
                    # continue

                ani_radius = float(sp_cluster_info['ani_radius'])

                results.append(run_sp_clustering_on_fastani_values(gid, gid_rep, df_fastani_subset, ani_radius, gid_suffix))

        log('Saving results')

        df = pd.DataFrame(results)
        self.save_hdf(df)
        return
