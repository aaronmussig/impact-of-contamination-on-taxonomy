import multiprocessing as mp
import os
import tempfile
from collections import defaultdict
from typing import Set

import pandas as pd
from chiton.fastani import fastani
from luigi import LocalTarget
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, DIR_OUT_FASTANI_INTER
from workflow.config import DIR_OUT_BATCH
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani.remove_gunc_failed_contigs_by_contamination_sp_cluster import \
    RemoveGuncFailedContigsByContaminationSpCluster
from workflow.fastani_interspecies.select_random_genomes_for_interspecies_ani_non_reps import \
    SelectRandomGenomesForInterspeciesAniNonReps
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty
import random

N_CPUS = mp.cpu_count()


class FastAniInterspeciesToGenusAniForRandomNonReps(LuigiTask):
    """
    Run FastANI on all sp cluster failed genomes to all species reps within the genus
    """

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'selected': SelectRandomGenomesForInterspeciesAniNonReps(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running FastANI interspecies genus ANI for random non reps', title=True)
        self.make_output_dirs()

        log('Loading selected gids')
        df_selected = self.input()['selected'].read() if not DEBUG else self.input()['selected'].read_cached()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        df_meta_reps = df_meta[df_meta['gtdb_representative'] == 't']

        log('Mapping genus to species representatives')
        d_genus_to_species_reps = defaultdict(set)
        for gid, row in tqdm(df_meta_reps.iterrows(), total=len(df_meta_reps)):
            d_genus_to_species_reps[row['genus']].add(gid)

        # Check if the batches have already been created
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        if DEBUG:
            batch_dir = '/tmp/b2a32x'
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir, exist_ok=True)

            log('Creating batch files')
            batch_paths = list()
            for genus, selected_gids in df_selected.itertuples():
                sp_reps = d_genus_to_species_reps[genus]
                cur_batch_path = os.path.join(batch_dir, f'batch_{genus}.tsv')
                with open(cur_batch_path, 'w') as f:
                    sp_reps_str = "|".join(sorted(sp_reps))
                    f.write(f'{selected_gids}\t{genus}\t{sp_reps_str}\n')
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


def fastani_worker(selected_gids: Set[str], genus, sp_reps, tmp_dir):

    path_results_out = os.path.join(DIR_OUT_FASTANI_INTER, 'for_random_non_reps', f'fastani_interspecies_to_{genus}.h5')
    # Stop early if this already exists
    if os.path.isfile(path_results_out):
        return

    # Generate the query gid paths
    d_qry_paths = dict()
    for qry_gid in selected_gids:
        d_qry_paths[qry_gid] = os.path.join(get_gid_root(qry_gid), f'{qry_gid}.fna')

    # Convert the rep set into their paths
    d_ref_paths = dict()
    d_path_to_gid = dict()
    for gid in sp_reps:
        d_ref_paths[gid] = os.path.join(tmp_dir, f'{gid}.fna')
        d_path_to_gid[d_ref_paths[gid]] = gid

    # Prepare for FastANI
    ani = fastani(query=list(d_qry_paths.values()),
                  reference=list(d_ref_paths.values()),
                  cpus=N_CPUS,
                  single_execution=False,
                  bidirectional=True,
                  show_progress=True)

    # Prepare the ANI results for output
    cur_results = list()
    d_ani = ani.as_dict()

    for qry_gid, qry_path in d_qry_paths.items():
        for ref_gid, ref_path in d_ref_paths.items():
            q_key = qry_path
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
            cur_results.append((qry_gid, ref_gid, ani, af))


    # Write the results to disk
    columns = ['query', 'reference', 'ani', 'af']
    path_results_tmp = os.path.join(tmp_dir, f'interspecies_ani_{genus}.h5')
    df = pd.DataFrame(cur_results, columns=columns)
    df.sort_values(by=['query', 'reference'], inplace=True)
    df.to_hdf(path_results_tmp, key='root', format='table', complevel=9, complib='blosc:lz4hc')
    move_file(path_results_tmp, path_results_out, checksum=True)
    return


def batch_worker(batch_path):
    log(f'Reading batch file: {batch_path}')
    cur_selected_gids, cur_genus, cur_sp_reps = None, None, None
    with open(batch_path) as f:
        for i, line in enumerate(f.readlines()):
            if i > 0:
                raise Exception('Malformed batch file')
            selected_gids, genus, sp_reps = line.strip().split('\t')
            sp_reps = frozenset(sp_reps.split('|'))
            cur_selected_gids = frozenset(selected_gids.split('|'))
            cur_genus = genus
            cur_sp_reps = sp_reps

    with tempfile.TemporaryDirectory() as tmp_dir:
        if DEBUG:
            tmp_dir = '/tmp/genomes'
            os.makedirs(tmp_dir, exist_ok=True)

        log(f'Caching genomes to scratch disk: {tmp_dir}')
        to_cache_queue = [(x, tmp_dir) for x in sp_reps]
        with mp.Pool(processes=min(mp.cpu_count(), 20)) as pool:
            list(tqdm(pool.imap_unordered(copy_genome_worker, to_cache_queue), total=len(to_cache_queue)))

        log(f'Running FastANI on {cur_genus}')
        fastani_worker(cur_selected_gids, cur_genus, cur_sp_reps, tmp_dir)
    return
