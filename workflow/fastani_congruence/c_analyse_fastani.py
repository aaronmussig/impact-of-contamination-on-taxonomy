import os

import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani_congruence.b_run_fastani import FastAniCongruenceRunFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastAniCongruenceAnalyseFastAni(LuigiTask):
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'ani': FastAniCongruenceRunFastAni(congruence=self.congruence, target_pct=self.target_pct,
                                               mash_k=self.mash_k, mash_s=self.mash_s),
            'meta': GtdbMetadataR207(),
            'sp_clusters': GtdbSpClustersR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                DIR_OUT_FASTANI_CONGRUENCE,
                f'results_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}.h5'
            ))

    def run(self):
        log(f'FastAniCongruenceAnalyseFastAni(p={self.target_pct}, c={self.congruence}, mash_k={self.mash_k}, mash_s={self.mash_s})',
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
            # print(gid)
            expected_species = df_meta.loc[gid, 'species']
            sp_cluster_info = df_sp_cluster.loc[expected_species]

            gid_rep = sp_cluster_info['rep_genome'][3:]

            try:
                df_fastani_subset = df_ani.loc[[gid]]
            except KeyError:
                # This happens if 0% of the genome was removed in the contig removal step
                results.append({
                    'gid': gid,
                    'new_sp_rep': None,
                    'ani': None,
                    'af': None,
                    'type': 'same_msa',
                    'same': True
                })
                continue

            ani_radius = float(sp_cluster_info['ani_radius'])

            results.append(run_sp_clustering_on_fastani_values(gid, gid_rep, df_fastani_subset, ani_radius))

        log('Saving results')

        df = pd.DataFrame(results)
        self.save_hdf(df)
        return
