import multiprocessing as mp
import os

import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTANI_INTER
from workflow.fastani.remove_gunc_failed_contigs_by_contamination_sp_cluster import \
    RemoveGuncFailedContigsByContaminationSpCluster
from workflow.fastani_interspecies.run_fastani_interspecies_to_genus_ani_for_sp_cluster_fail import \
    RunFastAniToGenusAniForSpClusterFail
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root

N_CPUS = mp.cpu_count()


class RunFastAniToGenusAniForSpClusterFailAgg(LuigiTask):
    """
    Run FastANI on all sp cluster failed genomes to all species reps within the genus
    """

    def requires(self):
        return {
            '_intergenus': RunFastAniToGenusAniForSpClusterFail(),
            'sp_cluster': RemoveGuncFailedContigsByContaminationSpCluster(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_FASTANI_INTER, 'fastani_fail_sp_clustering_interspecies_agg.h5'))

    def run(self):
        log('Running FastANI interspecies genus ANI for those which failed sp clustering (agg)', title=True)
        self.make_output_dirs()

        log('Loading species clustering failed genomes')
        df_sp = self.input()['sp_cluster'].read() if not DEBUG else self.input()['sp_cluster'].read_cached()
        fail_gids = frozenset({x[0] for x in df_sp[df_sp['same'] == False].index})
        log(f'Found {len(fail_gids):,}')

        queue = sorted(fail_gids)
        log(f'Processing queue: {len(queue):,}')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        df = pd.concat(results, ignore_index=True)
        self.save_hdf(df)


def worker(gid):
    gid_root = get_gid_root(gid)
    path_results_out = os.path.join(gid_root, 'fastani_interspecies_to_genus.h5')

    if not os.path.isfile(path_results_out):
        raise Exception(f'Not found: {gid}')

    df = pd.read_hdf(path_results_out)

    df['query'] = gid
    return df
