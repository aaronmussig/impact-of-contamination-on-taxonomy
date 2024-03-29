import multiprocessing as mp
import os
import tempfile

import pandas as pd
from Bio import SeqIO
from chiton.fastani import fastani
from magna.gunc import read_contig_assignments_tsv
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL
from workflow.config import DIR_OUT_BATCH
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty

N_CPUS = mp.cpu_count() // 2


N_CPUS = 15

class RemoveGuncFailedContigsByContamination(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'closest_reps': Closest100GenomesToRepresentative(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Removing contigs ordered by GUNC contamination (for failed only)', title=True)
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
            batch_dir = '/tmp/batch3'
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir, exist_ok=True)

            log('Creating queue')
            queue = list()
            for gid, row in tqdm(df_merged.iterrows(), total=len(df_merged)):
                cur_rep = row['gtdb_genome_representative'][3:]
                queue.append((
                    gid,
                    cur_rep,
                    d_rep_to_closest_gids[cur_rep],
                    row['source'],
                    row['taxonomic_level'],
                    row['domain'],
                ))
            queue = sorted(queue, key=lambda x: x[1])

            log('Creating batch files')
            batch_paths = list()
            for i, cur_batch in enumerate(iter_batches(queue, n=50)):
                cur_batch_path = os.path.join(batch_dir, f'batch_{i}.tsv')
                with open(cur_batch_path, 'w') as f:
                    for gid, cur_rep, closest_reps, gunc_source, gunc_taxlevel, domain in cur_batch:
                        f.write(f'{gid}\t{cur_rep}\t{closest_reps}\t'
                                f'{gunc_source}\t{gunc_taxlevel}\t{domain}\n')
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


def fastani_worker(gid, closest_rep_set, tmp_dir, gunc_source, max_css, domain):
    # Read the FASTA file
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')
    path_results_out = os.path.join(gid_root, 'fastani_gunc_failed_by_contamination.h5')

    # Stop early if this already exists
    if os.path.isfile(path_results_out):
        return

    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    if gunc_source == 'gtdb':
        path_contig_assign = os.path.join(gid_root, 'gunc_r95/gunc_output', f'{gid}.contig_assignments.tsv')
        gunc_domain = domain
    elif gunc_source == 'progenomes':
        path_contig_assign = os.path.join(gid_root, 'gunc_pro/gunc_output', f'{gid}.contig_assignments.tsv')
        if domain == 'd__Bacteria':
            gunc_domain = '2 Bacteria'
        elif domain == 'd__Archaea':
            gunc_domain = '2157 Archaea'
        else:
            raise ValueError(f'Unknown domain: {domain}')
    else:
        raise Exception(f'Unknown gunc source: {gunc_source}')
    df_contig_assign = read_contig_assignments_tsv(path_contig_assign)

    # Determine the percentage values at which this genome can have contigs removed
    taxon, tax_level = get_taxonomy_by_majority_vote_gunc(df_contig_assign, max_css, gunc_domain)
    d_pct_to_contigs_to_remove = contigs_to_remove_from_gunc(d_fna, df_contig_assign, taxon, tax_level)

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
    path_results_tmp = os.path.join(tmp_dir, 'fastani_gunc_failed_by_contamination.h5')
    df = pd.DataFrame(cur_results, columns=columns)
    df.to_hdf(path_results_tmp, key='root', format='table', complevel=9, complib='blosc:lz4hc')
    move_file(path_results_tmp, path_results_out, checksum=True)
    return


def batch_worker(batch_path):
    log(f'Reading batch file: {batch_path}')
    queue = list()
    with open(batch_path) as f:
        for line in f.readlines():
            gid, rep_gid, closest_reps, gunc_source, gunc_taxlevel, domain = line.strip().split('\t')
            closest_rep_set = frozenset(closest_reps.split('|'))
            if len(closest_rep_set) != 100:
                raise Exception(f'Expected 100 closest representatives, found {len(closest_rep_set)}')
            queue.append((gid, rep_gid, closest_rep_set, gunc_source, gunc_taxlevel, domain))
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
        for i, (gid, rep_gid, closest_rep_set, gunc_source, gunc_taxlevel, domain) in enumerate(queue):
            log(f'Processing {gid} {i}/{len(queue):,}')
            fastani_worker(gid, closest_rep_set, tmp_dir, gunc_source, gunc_taxlevel, domain)
    return
