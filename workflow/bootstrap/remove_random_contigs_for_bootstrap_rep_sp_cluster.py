import multiprocessing as mp
import os
import tempfile

import pandas as pd
from magna.util.disk import move_file
from tqdm import tqdm

from workflow.bootstrap.remove_random_contigs_for_bootstrap_rep import RemoveRandomContigsForBootstrapRep
from workflow.bootstrap.select_genomes_for_bootstrapping import SelectGenomesForBootstrapping
from workflow.config import DEBUG, DIR_OUT_BOOTSTRAP
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


def worker(job):
    gid, repeat, gid_i, gid_rep, gid_species, ani_radius = job

    path_out = os.path.join(DIR_OUT_BOOTSTRAP, 'random_contigs_sp_clustering', str(repeat), f'{gid}_{gid_i}.h5')
    if os.path.isfile(path_out):
        os.remove(path_out)
        # return pd.read_hdf(path_out)

    path_fastani = os.path.join(DIR_OUT_BOOTSTRAP, 'random_contigs', str(repeat), f'{gid}_{gid_i}.h5')
    if not os.path.isfile(path_fastani):
        raise Exception(f'{path_fastani} does not exist ({gid} {gid_i} {repeat})')
    try:
        df_fastani = pd.read_hdf(path_fastani)
    except ValueError as e:
        if os.path.getsize(path_fastani) != 1024:
            raise Exception(f'Error parsing: {path_fastani}')
        # This is an empty file, likely no contigs were removed - skip it
        return None

    pct_values = frozenset(int(x) for x in df_fastani['pct'].unique())

    out = list()
    for pct_value in sorted(pct_values):
        df_fastani_subset = df_fastani[df_fastani['pct'] == pct_value]
        result = run_sp_clustering_on_fastani_values(gid, gid_rep, df_fastani_subset, ani_radius)
        result['pct'] = pct_value
        result['repeat'] = repeat
        result['gid_i'] = gid_i
        out.append(result)

    # Save intermediate file to disk
    os.makedirs(os.path.dirname(path_out), exist_ok=True)
    df = pd.DataFrame(out)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, 'out.h5')
        df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')
        move_file(tmp_path, path_out, checksum=True)
    return df


class RemoveRandomContigsForBootstrapRepSpCluster(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            '_remove_contigs': RemoveRandomContigsForBootstrapRep(),
            'sp_clusters': GtdbSpClustersR207(),
            'meta': GtdbMetadataR207(),
            'bootstrap_gids': SelectGenomesForBootstrapping(),
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(DIR_OUT_BOOTSTRAP, 'remove_random_contigs_for_bootstrap_rep_sp_cluster.h5'))

    def run(self):
        log('Running species clustering on gunc failed contig by bootstrap random', title=True)
        self.make_output_dirs()

        log('Loading bootstrap genomes')
        df_bootstrap = self.input()['bootstrap_gids'].read() if not DEBUG else self.input()[
            'bootstrap_gids'].read_cached()

        log('Loading GTDB metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading GTDB species clusters')
        df_sp_clusters = self.input()['sp_clusters'].read() if not DEBUG else self.input()['sp_clusters'].read_cached()

        log('Creating queue')
        queue = list()
        for repeat, row in tqdm(df_bootstrap.iterrows(), total=len(df_bootstrap)):
            gids = row['genomes'].split('|')
            for gid_i, gid in enumerate(gids):
                gid_rep = df_meta.loc[gid, 'gtdb_genome_representative'][3:]
                gid_species = df_meta.loc[gid, 'species']
                ani_radius = float(df_sp_clusters.loc[gid_species, 'ani_radius'])
                queue.append((
                    gid,
                    repeat,
                    gid_i,
                    gid_rep,
                    gid_species,
                    ani_radius
                ))
        log(f'Queue size: {len(queue):,}')

        log('Running species clustering')
        if DEBUG:
            results = [worker(x) for x in tqdm(queue)]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log(f'Concatenating dataframe: {len(results):,}')
        results = [x for x in results if x is not None]
        log(f'Results size: {len(results):,}')
        df = pd.concat(results, ignore_index=True)

        log('Saving dataframe')
        if not DEBUG:
            self.save_hdf(df)
