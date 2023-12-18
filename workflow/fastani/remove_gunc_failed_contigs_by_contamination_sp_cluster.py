import multiprocessing as mp
import os
import tempfile

import pandas as pd
from magna.util.disk import move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTANI
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani.remove_gunc_failed_contigs_by_contamination_aggregate import \
    RemoveGuncFailedContigsByContaminationAggregate
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def worker(job):
    gid, df_sp_clusters, df_meta = job
    gid_root = get_gid_root(gid)

    path_out = os.path.join(gid_root, 'fastani_gunc_failed_by_contamination_sp_clustering.h5')
    if os.path.isfile(path_out):
        return pd.read_hdf(path_out)

    df_fastani = pd.read_hdf(os.path.join(gid_root, 'fastani_gunc_failed_by_contamination.h5'))

    gid_rep = df_meta.loc[gid, 'gtdb_genome_representative'][3:]
    gid_species = df_meta.loc[gid, 'species']

    pct_values = frozenset(int(x) for x in df_fastani['pct'].unique())

    ani_radius = float(df_sp_clusters.loc[gid_species, 'ani_radius'])

    out = list()
    for pct_value in pct_values:
        df_fastani_subset = df_fastani[df_fastani['pct'] == pct_value]
        result = run_sp_clustering_on_fastani_values(gid, gid_rep, df_fastani_subset, ani_radius)
        result['pct'] = pct_value
        out.append(result)

    # Save intermediate file to disk
    df = pd.DataFrame(out)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, 'out.h5')
        df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')
        move_file(tmp_path, path_out, checksum=True)
    return df


class RemoveGuncFailedContigsByContaminationSpCluster(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'fastani': RemoveGuncFailedContigsByContaminationAggregate(),
            'sp_clusters': GtdbSpClustersR207(),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(DIR_OUT_FASTANI, 'remove_gunc_failed_contigs_by_contamination_sp_clustering.h5'))

    def run(self):
        log('Running species clustering on gunc failed contig by contamination', title=True)
        self.make_output_dirs()

        log('Loading GTDB metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Loading GTDB species clusters')
        df_sp_clusters = self.input()['sp_clusters'].read() if not DEBUG else self.input()['sp_clusters'].read_cached()

        log('Creating queue')
        queue = list()
        for gid in tqdm(sorted(df_css.index), total=len(df_css.index)):
            queue.append((gid, df_sp_clusters, df_meta))
            if DEBUG and len(queue) > 10:
                break
        log(f'Queue size: {len(queue):,}')

        log('Running species clustering')
        results = [worker(x) for x in tqdm(queue)]
        # if DEBUG:
        #     results = [worker(x) for x in tqdm(queue)]
        # else:
        #     with mp.Pool(processes=10) as pool:
        #         results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log('Concatenating dataframe')
        df = pd.concat(results, ignore_index=True)

        log('Saving dataframe')
        if not DEBUG:
            self.save_hdf(df, index=['gid', 'pct'])
