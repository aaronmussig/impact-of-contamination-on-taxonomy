import multiprocessing as mp

import pandas as pd

from workflow.config import DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani.run_fastani_on_interspecies_aggregate import RunFastAniOnInterspeciesAniAggregate
from workflow.fastani_interspecies.run_fastani_interspecies_to_genus_ani_for_random_non_reps_agg import \
    FastAniInterspeciesToGenusAniForRandomNonRepsAgg
from workflow.fastani_interspecies.run_fastani_interspecies_to_genus_ani_for_sp_cluster_fail_agg import \
    RunFastAniToGenusAniForSpClusterFailAgg
from workflow.model.luigi import LuigiTask
from workflow.util.log import log

N_CPUS = mp.cpu_count()

N_RANDOM = 10
import matplotlib.pyplot as plt
import seaborn as sns


class FastAniInterspeciesPlots(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'fail': RunFastAniToGenusAniForSpClusterFailAgg(),
            'inter_non': FastAniInterspeciesToGenusAniForRandomNonRepsAgg(),
            'inter_rep': RunFastAniOnInterspeciesAniAggregate(),
        }

    def run(self):
        log('Creating plots for FastANI interspecies', title=True)

        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        s_rep_gids = frozenset(df_meta[df_meta['gtdb_representative'] == 't'].index)

        # Load the gids that failed during species clustering, and get their distances to reps
        log('Loading fail')
        df_fail = self.input()['fail'].read() if not DEBUG else self.input()['fail'].read_cached()
        df_fail = df_fail[(df_fail['af'] >= 0.5) & (df_fail['query'] != df_fail['reference'])]
        log(df_fail.shape)

        # Load the species rep inter-genus calculation
        log('Loading reps')
        df_inter_rep = self.input()['inter_rep'].read() if not DEBUG else self.input()['inter_rep'].read_cached()
        df_inter_rep = df_inter_rep[(df_inter_rep['af'] >= 0.5) & (df_inter_rep['query'] != df_inter_rep['reference'])]
        log(df_inter_rep.shape)

        # Load the species non rep inter-genus calc
        log('Loading non-reps')
        df_inter_non = self.input()['inter_non'].read() if not DEBUG else self.input()['inter_non'].read_cached()
        df_inter_non = df_inter_non[(df_inter_non['af'] >= 0.5) & (df_inter_non['query'] != df_inter_non['reference'])]
        log(df_inter_non.shape)

        rows = list()
        for row in df_fail.itertuples():
            rows.append({
                'ani': row.ani,
                'is_rep': row.query in s_rep_gids,
                'source': 'fail'
            })
        for row in df_inter_rep.itertuples():
            rows.append({
                'ani': row.ani,
                'is_rep': row.query in s_rep_gids,
                'source': 'all',
            })
        for row in df_inter_non.itertuples():
            rows.append({
                'ani': row.ani,
                'is_rep': row.query in s_rep_gids,
                'source': 'all',
            })
        df = pd.DataFrame(rows)

        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
        # sns.histplot(data=df[df['source'] == 'fail'], x='ani', element='poly',
        #              stat='percent', common_norm=True, ax=ax1)
        # sns.histplot(data=df[df['source'] == 'all'], x='ani',  element='poly',
        #              stat='percent', common_norm=True, ax=ax2)
        # plt.show()


        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
        # plt.rcParams.update({'font.size': 13})
        # plt.rcParams['svg.fonttype'] = 'none'
        #
        # sns.histplot(data=df[df['source'] == 'fail'], x='ani',
        #              bins=list(range(75, 100)),
        #              stat='percent', common_norm=True, ax=ax2)
        # sns.histplot(data=df[df['source'] == 'all'], x='ani',
        #              bins=list(range(75, 100)),
        #              stat='percent', common_norm=True, ax=ax1)
        # plt.show()

        fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))
        plt.rcParams.update({'font.size': 13})
        plt.rcParams['svg.fonttype'] = 'none'

        sns.histplot(data=df, x='ani',
                     bins=list(range(75, 100)), hue='source', multiple='stack',
                     stat='percent', common_norm=False, ax=ax1)

        plt.savefig('/tmp/new_ani.svg')
        plt.show()

        return
