import multiprocessing as mp
import os
import tempfile

import pandas as pd
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL
from workflow.config import DIR_OUT_BATCH
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.gunc.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty

N_CPUS = mp.cpu_count() // 2


class RunFastAniOnGuncFailedBaseCase(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'closest_reps': Closest100GenomesToRepresentative(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running FastANI (base case=0) on GUNC failed genomes', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        log('Loading closest representatives dataframe')
        df_closest_reps = self.input()['closest_reps'].read() if not DEBUG else self.input()[
            'closest_reps'].read_cached()
        d_rep_to_closest_gids = df_closest_reps['closest_representatives'].to_dict()
        log(f'Found {len(d_rep_to_closest_gids)} representatives')

        # Check if the batches have already been created
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        if DEBUG:
            batch_dir = '/tmp/batc4h3'
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir, exist_ok=True)

            log('Creating queue')
            queue = list()
            for gid, row in tqdm(df_merged.iterrows(), total=len(df_merged)):
                cur_rep = row['gtdb_genome_representative'][3:]
                queue.append((
                    gid,
                    cur_rep,
                    d_rep_to_closest_gids[cur_rep]
                ))
            queue = sorted(queue, key=lambda x: x[1])

            log('Creating batch files')
            batch_paths = list()
            for i, cur_batch in enumerate(iter_batches(queue, n=50 if not DEBUG else 2)):
                cur_batch_path = os.path.join(batch_dir, f'batch_{i}.tsv')
                with open(cur_batch_path, 'w') as f:
                    for gid, cur_rep, closest_reps in cur_batch:
                        f.write(f'{gid}\t{cur_rep}\t{closest_reps}\n')
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


def fastani_worker(gid, closest_rep_set, tmp_dir):
    # Read the FASTA file
    gid_root = get_gid_root(gid)
    qry_path = os.path.join(gid_root, f'{gid}.fna')
    path_results_out = os.path.join(gid_root, 'fastani_gunc_failed_pct_0.h5')

    # Stop early if this already exists
    if os.path.isfile(path_results_out):
        return

    # Convert the rep set into their paths
    d_ref_paths = {x: os.path.join(tmp_dir, f'{x}.fna') for x in closest_rep_set}

    # Prepare for FastANI
    ani = fastani(query=qry_path,
                  reference=list(d_ref_paths.values()),
                  cpus=N_CPUS,
                  single_execution=False,
                  bidirectional=True,
                  show_progress=False)

    # Prepare the ANI results for output
    cur_results = list()
    d_ani = ani.as_dict()
    for ref_gid in d_ref_paths.keys():
        q_key = qry_path
        r_key = d_ref_paths[ref_gid]
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
        cur_results.append((ref_gid, 0, ani, af))

    # Write the results to disk
    columns = ['reference', 'pct', 'ani', 'af']
    path_results_tmp = os.path.join(tmp_dir, 'fastani_gunc_failed_pct_0.h5')
    df = pd.DataFrame(cur_results, columns=columns)
    df.to_hdf(path_results_tmp, key='root', format='table', complevel=9, complib='blosc:lz4hc')
    move_file(path_results_tmp, path_results_out, checksum=True)
    return


def batch_worker(batch_path):
    log(f'Reading batch file: {batch_path}')
    queue = list()
    with open(batch_path) as f:
        for line in f.readlines():
            gid, rep_gid, closest_reps = line.strip().split('\t')
            closest_rep_set = frozenset(closest_reps.split('|'))
            if len(closest_rep_set) != 100:
                raise Exception(f'Expected 100 closest representatives, found {len(closest_rep_set)}')
            queue.append((gid, rep_gid, closest_rep_set))
    log(f'Found {len(queue):,} entries in batch file')

    # if DEBUG:
    #     queue = [x for x in queue if x[0] == 'GCA_001404855.1']
    #     if len(queue) == 0:
    #         return

    log('Finding unique representative gids for this batch')
    unq_reps = set()
    [unq_reps.update(x[2]) for x in queue]
    unq_reps = sorted(unq_reps)
    log(f'Found {len(unq_reps):,} unique representatives')

    with tempfile.TemporaryDirectory() as tmp_dir:
        if DEBUG:
            tmp_dir = '/tmp/genomes'
            os.makedirs(tmp_dir, exist_ok=True)

        log(f'Caching genomes to scratch disk: {tmp_dir}')
        to_cache_queue = [(x, tmp_dir) for x in unq_reps]
        with mp.Pool(processes=min(mp.cpu_count(), 20)) as pool:
            list(tqdm(pool.imap_unordered(copy_genome_worker, to_cache_queue), total=len(to_cache_queue)))

        log(f'Running FastANI on {len(queue):,} genomes')
        for i, (gid, rep_gid, closest_rep_set) in enumerate(queue):
            log(f'Processing {gid} {i}/{len(queue):,}')
            fastani_worker(gid, closest_rep_set, tmp_dir)
    return
