import multiprocessing as mp
import os
import tempfile
from collections import defaultdict

import pandas as pd
from chiton.fastani import fastani
from luigi import LocalTarget
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL
from workflow.config import DIR_OUT_BATCH
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani.remove_gunc_failed_contigs_by_contamination_sp_cluster import \
    RemoveGuncFailedContigsByContaminationSpCluster
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty

N_CPUS = mp.cpu_count()


class RunFastAniToGenusAniForSpClusterFail(LuigiTask):
    """
    Run FastANI on all sp cluster failed genomes to all species reps within the genus
    """

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'sp_cluster': RemoveGuncFailedContigsByContaminationSpCluster(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running FastANI interspecies genus ANI for those which failed sp clustering', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        df_meta_reps = df_meta[df_meta['gtdb_representative'] == 't']

        log('Loading species clustering failed genomes')
        df_sp = self.input()['sp_cluster'].read() if not DEBUG else self.input()['sp_cluster'].read_cached()
        fail_gids = frozenset({x[0] for x in df_sp[df_sp['same'] == False].index})
        log(f'Found {len(fail_gids):,}')

        d_fail_genus_to_gids = defaultdict(set)
        for fail_gid in fail_gids:
            d_fail_genus_to_gids[df_meta.loc[fail_gid, 'genus']].add(fail_gid)
        log(f'Belonging to {len(d_fail_genus_to_gids):,} genera')

        log('Mapping genus to species representatives')
        d_genus_to_species_reps = defaultdict(set)
        for gid, row in tqdm(df_meta_reps.iterrows(), total=len(df_meta_reps)):
            d_genus_to_species_reps[row['genus']].add(gid)

        # Check if the batches have already been created
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        if DEBUG:
            batch_dir = '/tmp/ba2223tc14h33'
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir, exist_ok=True)

            log('Creating queue')
            queue = list()
            for genus, fail_gids in d_fail_genus_to_gids.items():
                sp_reps = d_genus_to_species_reps[genus]
                for fail_gid in fail_gids:
                    queue.append((fail_gid, genus, sorted(sp_reps)))

            log('Creating batch files')
            batch_paths = list()
            for fail_gid, genus, sp_reps in sorted(queue, key=lambda x: x[0]):
                cur_batch_path = os.path.join(batch_dir, f'batch_{fail_gid}.tsv')
                with open(cur_batch_path, 'w') as f:
                    f.write(f'{fail_gid}\t{genus}\t{"|".join(sp_reps)}\n')
                batch_paths.append((cur_batch_path,))

            # Run each of the batches
            if not DEBUG:
                log('Submitting batches to RQ...')
                submit_jobs_to_rq(fn=batch_worker, q_args=batch_paths, queue_name=self.fqn)
        else:
            batch_paths = [(os.path.join(batch_dir, x),) for x in os.listdir(batch_dir)]

        if DEBUG:
            log('Starting workers (single-threaded)')
            [batch_worker(x[0]) for x in batch_paths]
        else:
            # Run each of the batches
            log('Waiting until RQ queue is empty...')
            rq_wait_for_queue_empty(self.fqn)
            self.write_sentinel()


def copy_genome_worker(job):
    gid, tmp_dir = job
    gid_path_srv = os.path.join(get_gid_root(gid), f'{gid}.fna')
    gid_path_tmp = os.path.join(tmp_dir, f'{gid}.fna')
    if not os.path.isfile(gid_path_tmp):
        copy_file(gid_path_srv, gid_path_tmp)


def fastani_worker(fail_gid, genus, sp_reps, tmp_dir):
    gid_root = get_gid_root(fail_gid)
    path_results_out = os.path.join(gid_root, 'fastani_interspecies_to_genus.h5')
    query_fna_path = os.path.join(gid_root, f'{fail_gid}.fna')

    # Stop early if this already exists
    if os.path.isfile(path_results_out):
        return

    # Convert the rep set into their paths
    d_ref_paths = dict()
    d_path_to_gid = dict()
    gids = list()
    for gid in sp_reps:
        d_ref_paths[gid] = os.path.join(tmp_dir, f'{gid}.fna')
        d_path_to_gid[d_ref_paths[gid]] = gid
        gids.append(gid)

    # Prepare for FastANI
    ani = fastani(query=query_fna_path,
                  reference=list(d_ref_paths.values()),
                  cpus=N_CPUS,
                  single_execution=False,
                  bidirectional=True,
                  show_progress=True)

    # Prepare the ANI results for output
    cur_results = list()
    d_ani = ani.as_dict()

    for ref_gid, ref_path in d_ref_paths.items():
        q_key = query_fna_path
        r_key = ref_path
        qvr = d_ani[q_key][r_key]
        rvq = d_ani[r_key][q_key]

        if qvr is not None and rvq is not None:
            ani = max(qvr.ani, rvq.ani)
            af = max(qvr.align_frac, rvq.align_frac)
        elif qvr is not None and rvq is None:
            ani = qvr.ani
            af = qvr.align_frac
        elif qvr is None and rvq is not None:
            ani = rvq.ani
            af = rvq.align_frac
        else:
            ani = 0
            af = 0
        cur_results.append((ref_gid, ani, af))

    # Write the results to disk
    columns = ['reference', 'ani', 'af']
    path_results_tmp = os.path.join(tmp_dir, f'interspecies_ani_{fail_gid}.h5')
    df = pd.DataFrame(cur_results, columns=columns)
    df.sort_values(by='reference', inplace=True)
    df.to_hdf(path_results_tmp, key='root', format='table', complevel=9, complib='blosc:lz4hc')
    move_file(path_results_tmp, path_results_out, checksum=True)
    return


def batch_worker(batch_path):
    log(f'Reading batch file: {batch_path}')
    queue = list()
    gids_to_cache = set()
    with open(batch_path) as f:
        for line in f.readlines():
            fail_gid, genus, sp_reps = line.strip().split('\t')
            sp_reps = frozenset(sp_reps.split('|'))
            gids_to_cache.update(sp_reps)
            queue.append((fail_gid, genus, sp_reps))
    log(f'Found {len(queue):,} entries in batch file')

    with tempfile.TemporaryDirectory() as tmp_dir:
        if DEBUG:
            tmp_dir = '/tmp/genomes'
            os.makedirs(tmp_dir, exist_ok=True)

        log(f'Caching genomes to scratch disk: {tmp_dir}')
        to_cache_queue = [(x, tmp_dir) for x in gids_to_cache]
        with mp.Pool(processes=min(mp.cpu_count(), 20)) as pool:
            list(tqdm(pool.imap_unordered(copy_genome_worker, to_cache_queue), total=len(to_cache_queue)))

        log(f'Running FastANI on {len(queue):,} genomes')
        for fail_gid, genus, sp_reps in queue:
            fastani_worker(fail_gid, genus, sp_reps, tmp_dir)
    return
