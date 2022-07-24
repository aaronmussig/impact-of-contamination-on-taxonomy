import multiprocessing as mp
import os
import tempfile
from typing import Dict, Set

import pandas as pd
from Bio import SeqIO
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_BACKGROUND, PCT_VALUES
from workflow.config import DIR_OUT_BATCH
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty

N_CPUS = mp.cpu_count()


def get_pct_to_remove_contigs_ordered_by_size(d_fna: Dict[str, str]) -> Dict[int, Set[str]]:
    d_contig_to_len = {k: len(v) for k, v in d_fna.items()}
    genome_len = sum(d_contig_to_len.values())
    d_contig_to_len_sorted = sorted(d_contig_to_len.items(), key=lambda x: (x[1], x[0]))
    d_pct_to_contigs = dict()
    for pct in PCT_VALUES:
        contigs_to_remove = set()
        total_len_removed = 0
        for contig, contig_len in d_contig_to_len_sorted:
            contig_len = d_contig_to_len[contig]
            new_len_removed = total_len_removed + contig_len
            new_pct_removed = new_len_removed / genome_len * 100

            # Don't add it as it will be over the threshold
            if new_pct_removed > pct:
                continue

            contigs_to_remove.add(contig)
            total_len_removed += contig_len

        d_pct_to_contigs[pct] = contigs_to_remove

    # De-duplicate the data set to only those percentages where a unique number of contigs were removed
    last_contigs_removed = set()
    d_pct_to_contigs_dedup = dict()
    for pct, contigs_removed in d_pct_to_contigs.items():
        if contigs_removed == last_contigs_removed:
            continue
        d_pct_to_contigs_dedup[pct] = contigs_removed
        last_contigs_removed = contigs_removed
    return d_pct_to_contigs_dedup


class BackgroundRemoveContigsOrderedBySize(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'closest_reps': Closest100GenomesToRepresentative()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_BACKGROUND, self.fqn))

    def run(self):
        log('Removing contigs ordered by length for background comparison', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading closest representatives dataframe')
        df_closest_reps = self.input()['closest_reps'].read() if not DEBUG else self.input()[
            'closest_reps'].read_cached()
        d_rep_to_closest_gids = df_closest_reps['closest_representatives'].to_dict()
        log(f'Found {len(d_rep_to_closest_gids)} representatives')

        # Check if the batches have already been created
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir, exist_ok=True)

            log('Creating queue')
            queue = list()
            for gid, row in tqdm(df_meta.iterrows(), total=len(df_meta)):
                cur_rep = row['gtdb_genome_representative'][3:]
                queue.append((gid, cur_rep, d_rep_to_closest_gids[cur_rep]))
                if DEBUG and len(queue) >= 1:
                    break
            queue = sorted(queue, key=lambda x: x[1])

            log('Creating batch files')
            batch_paths = list()
            for i, cur_batch in enumerate(iter_batches(queue, n=50)):
                cur_batch_path = os.path.join(batch_dir, f'batch_{i}.tsv')
                with open(cur_batch_path, 'w') as f:
                    for gid, cur_rep, closest_reps in cur_batch:
                        f.write(f'{gid}\t{cur_rep}\t{closest_reps}\n')
                batch_paths.append((cur_batch_path,))

            # Run each of the batches
            log('Submitting batches to RQ...')
            submit_jobs_to_rq(fn=batch_worker, q_args=batch_paths, queue_name=self.fqn)
        else:
            batch_paths = [(os.path.join(batch_dir, x), ) for x in os.listdir(batch_dir)]

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
    gid_fna = os.path.join(gid_root, f'{gid}.fna')
    path_results_out = os.path.join(gid_root, 'fastani_background_by_contig_size.h5')

    # Stop early if this already exists
    if os.path.isfile(path_results_out):
        return

    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    # Determine the percentage values at which this genome can have contigs removed
    d_pct_to_contigs_to_remove = get_pct_to_remove_contigs_ordered_by_size(d_fna)

    log(f'Running FastANI on {len(d_pct_to_contigs_to_remove):,} pct values')

    # Convert the rep set into their paths
    d_ref_paths = {x: os.path.join(tmp_dir, f'{x}.fna') for x in closest_rep_set}

    # Iterate over each percentage change
    cur_results = list()
    for pct, contigs_to_remove in tqdm(d_pct_to_contigs_to_remove.items()):

        # Write the fasta file without those contigs to disk
        qry_path = os.path.join(tmp_dir, f'{gid}_{pct}.fna')
        with open(qry_path, 'w') as f:
            SeqIO.write([x[1] for x in d_fna.items() if x[0] not in contigs_to_remove], f, 'fasta')

        # Prepare for FastANI
        ani = fastani(query=qry_path,
                      reference=list(d_ref_paths.values()),
                      cpus=N_CPUS,
                      single_execution=False,
                      bidirectional=True,
                      show_progress=False,
                      tmp_root='/dev/shm' if not DEBUG else tempfile.gettempdir())

        # Prepare the ANI results for output
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
            cur_results.append((ref_gid, pct, ani, af))

    # Write the results to disk
    columns = ['reference', 'pct', 'ani', 'af']
    path_results_tmp = os.path.join(tmp_dir, 'fastani_background_by_contig_size.h5')
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
