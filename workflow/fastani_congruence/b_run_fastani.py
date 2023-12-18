import multiprocessing as mp
import os
import tempfile
from collections import defaultdict
from typing import Set

import luigi
import pandas as pd
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE, DEBUG, DIR_OUT_BATCH, DIR_CACHE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani_congruence.a_run_mash_on_failed import FastAniCongruenceRunMash
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.gunc_helper.genome_pct_congruence_contigs_removed import GenomePctCongruenceContigsRemoved
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_wait_for_queue_empty, rq_and_wait


class FastAniCongruenceRunFastAni(LuigiTask):
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'mash': FastAniCongruenceRunMash(
                target_pct=self.target_pct,
                congruence=self.congruence,
                mash_k=self.mash_k,
                mash_s=self.mash_s
            ),
            'contigs_removed': GenomePctCongruenceContigsRemoved(congruence=self.congruence, target_pct=self.target_pct),
            'meta': GtdbMetadataR207(),
            'closest_reps': Closest100GenomesToRepresentative()
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                DIR_OUT_FASTANI_CONGRUENCE,
                f'ani_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}.h5'
            ))

    def run(self):
        log(f'FastAniCongruenceRunFastAni(p={self.target_pct}, c={self.congruence})', title=True)
        self.make_output_dirs()

        batch_size = 25

        batch_root = os.path.join(DIR_OUT_BATCH,
                                  f'{self.fqn}_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}')
        if not os.path.isdir(batch_root) or True:
            log('Creating batch paths')
            os.makedirs(batch_root, exist_ok=True)
            batch_paths = self.create_batch_files(batch_root, batch_size)
            rq_submit = True
        else:
            log('Loading batch paths')
            batch_paths = self.load_batch_files(batch_root)
            rq_submit = False

        if rq_submit:
            queue = [(x, self.target_pct, self.congruence, self.mash_k, self.mash_s) for x in batch_paths]
            log('Submitting batches to RQ...')
            job_id = f'{self.fqn}_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}'
            log(f'Job ID: {job_id}')
            rq_and_wait(job_id=job_id, fn=run_on_batch_file, q_args=queue, queue_name=job_id)
            self.write_sentinel()

            # Run each of the batches
            log('Waiting until RQ queue is empty...')
            rq_wait_for_queue_empty(job_id)

        # Collecting the results
        log('Collecting results')


        log('Reading batch files')
        queue = list()
        for job, _ in [read_batch_file(x) for x in batch_paths]:
            for qry_gid in job.keys():
                queue.append((qry_gid, self.congruence, self.target_pct, self.mash_k, self.mash_s))

            if DEBUG:
                break


        log('Raeding dataframes')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            dfs = list(tqdm(pool.imap_unordered(collect_result_worker, queue), total=len(queue)))

        log('Concat dfs')
        df = pd.concat(dfs, ignore_index=True)
        df.sort_values(by=['query', 'ani', 'af', 'ref'], ascending=[True, False, False, True], inplace=True, ignore_index=True)
        log('Writing output')
        if not DEBUG:
            self.save_hdf(df)

        return

    def create_batch_files(self, batch_root, batch_size):


        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        df_meta_rep = df_meta[df_meta['gtdb_representative'] == 't']
        d_rank_to_rep_gids = defaultdict(set)
        for cur_rank in ('phylum', 'class', 'order', 'family', 'genus'):
            for gid, taxon in df_meta_rep[cur_rank].to_dict().items():
                d_rank_to_rep_gids[taxon].add(gid)

        log('Loading contigs removed from failed genomes')
        df_removed = self.input()['contigs_removed'].maybe_read_cached()
        d_gid_to_contigs_removed = dict()
        for row in df_removed.itertuples():
            d_gid_to_contigs_removed[row.Index] = row.contigs_removed.split('|')

        log('Allocating genomes into batches')
        d_batch_id_to_gids = defaultdict(set)
        for batch_id, batch_gids in enumerate(iter_batches(sorted(d_gid_to_contigs_removed.keys()), n=batch_size)):
            for batch_gid in batch_gids:
                d_batch_id_to_gids[batch_id].add(batch_gid)

        d_qry_to_refs = defaultdict(set)
        log('Adding genomes from the expected taxa')
        for qry_gid in tqdm(d_gid_to_contigs_removed.keys()):
            cur_meta_row = df_meta.loc[qry_gid]
            for cur_rank in ('genus', 'family', 'order', 'class', 'phylum'):
                cur_taxon = cur_meta_row[cur_rank]
                gids_in_taxon = d_rank_to_rep_gids[cur_taxon]
                if cur_rank != 'genus':
                    if len(gids_in_taxon) > 100:
                        break
                    else:
                        d_qry_to_refs[qry_gid].update(gids_in_taxon)
                else:
                    d_qry_to_refs[qry_gid].update(gids_in_taxon)

        log('Adding closest representative genomes (by patristic distance)')
        df_closest_reps = self.input()['closest_reps'].maybe_read_cached()
        for qry_gid in tqdm(d_gid_to_contigs_removed.keys()):
            cur_sp_rep_gid = df_meta.loc[qry_gid, 'gtdb_genome_representative'][3:]
            closest_reps = frozenset(df_closest_reps.loc[cur_sp_rep_gid, 'closest_representatives'].split('|'))
            d_qry_to_refs[qry_gid].update(closest_reps)

        log('Loading Mash scores')
        df_mash = self.input()['mash'].maybe_read_cached()
        for row in tqdm(df_mash.itertuples(), total=len(df_mash)):
            d_qry_to_refs[row.query].add(row.ref)

        log('Creating batch files')
        os.makedirs(batch_root, exist_ok=True)
        batch_paths = list()
        for batch_id, batch_gids in d_batch_id_to_gids.items():
            batch_path = os.path.join(batch_root, f'{str(batch_id)}.tsv')
            batch_paths.append(batch_path)
            with open(batch_path, 'w') as f:
                for batch_gid in batch_gids:
                    rep_gids_str = '|'.join(d_qry_to_refs[batch_gid])
                    if len(rep_gids_str) < 2:
                        print(f'{batch_gid} has no mash hits')
                    else:
                        contigs_to_remove_str = '|'.join(d_gid_to_contigs_removed[batch_gid])
                        if len(contigs_to_remove_str) < 2:
                            print(f'{batch_gid} has no contigs to remove')
                        else:
                            f.write(f'{batch_gid}\t{rep_gids_str}\t{contigs_to_remove_str}\n')
        return batch_paths

    def load_batch_files(self, batch_root):
        out = list()
        for file in os.listdir(batch_root):
            if file.endswith('.tsv'):
                out.append(os.path.join(batch_root, file))
        return out


def run_fastani_on_jobs(qry_gid: str, ref_gids: Set[str], tmp_dir, tmp_ref_dir, path_out: str):
    # if os.path.exists(path_out):
    #     if DEBUG:
    #         df = pd.read_hdf(path_out)
    #         print()
    #     return

    path_fna_tmp = os.path.join(get_gid_root(qry_gid, root=tmp_dir), f'{qry_gid}.fna')
    ref_gid_paths = {x: os.path.join(get_gid_root(x, root=tmp_ref_dir), f'{x}.fna') for x in ref_gids}

    ani = fastani(query=path_fna_tmp,
                  reference=list(ref_gid_paths.values()),
                  cpus=mp.cpu_count(),
                  single_execution=False,
                  bidirectional=True,
                  show_progress=False)

    # Prepare the ANI results for output
    d_ani = ani.as_dict()
    out = list()
    for ref_gid in ref_gid_paths.keys():
        q_key = path_fna_tmp
        r_key = ref_gid_paths[ref_gid]
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

        out.append({
            'query': qry_gid,
            'ref': ref_gid,
            'ani': ani,
            'af': af,
            'ani_qvr': qvr.ani if qvr is not None else 0,
            'ani_rvq': rvq.ani if rvq is not None else 0,
            'af_qvr': qvr.align_frac if qvr is not None else 0,
            'af_rvq': rvq.align_frac if rvq is not None else 0,
        })

    df = pd.DataFrame(out)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, 'out.h5')
        df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')
        move_file(tmp_path, path_out, checksum=True)

    return


def write_genome_fna_worker(job):
    gid, contigs_removed, tmp_dir = job

    genome = Genome(gid)

    path_fna_out = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}.fna')
    os.makedirs(os.path.dirname(path_fna_out), exist_ok=True)

    # Write the modified FNA file to disk
    with open(path_fna_out, 'w') as f:
        for contig, seq in genome.d_fna.items():
            if contig not in contigs_removed:
                f.write(f'>{contig}\n{seq}\n')

    return


def genome_cache_worker(job):
    gid, tmp_dir = job
    org_gid = os.path.join(get_gid_root(gid), f'{gid}.fna')
    new_path = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}.fna')

    os.makedirs(os.path.dirname(new_path), exist_ok=True)
    if not os.path.exists(new_path):
        copy_file(org_gid, new_path)
    return


def read_batch_file(path):
    ref_gids = set()
    d_gid_to_jobs = dict()
    with open(path, 'r') as f:
        for line in f.readlines():
            gid, rep_gids_str, contigs_removed_str = line.strip().split('\t')
            rep_gids = frozenset(rep_gids_str.split('|'))
            contigs_removed = frozenset(contigs_removed_str.split('|'))
            d_gid_to_jobs[gid] = (rep_gids, contigs_removed)
            ref_gids.update(rep_gids)
    return d_gid_to_jobs, ref_gids


def run_on_batch_file(batch_path, target_pct, congruence, mash_k, mash_s):
    log(f'Processing batch file: {batch_path}')
    d_gid_to_jobs, ref_gids = read_batch_file(batch_path)

    rep_gid_dir = os.path.join(DIR_CACHE, 'genomes')

    with tempfile.TemporaryDirectory() as tmp_dir:

        log('Getting the pruned genome file for all failed')
        queue = list()
        queued_gids = set()
        for gid, (rep_gids, contigs_removed) in tqdm(d_gid_to_jobs.items()):
            queue.append((
                gid,
                contigs_removed,
                tmp_dir,
            ))
            queued_gids.add(gid)

        log('Processing queue')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            list(tqdm(pool.imap_unordered(write_genome_fna_worker, queue), total=len(queue), smoothing=0.05))

        log('Generating queue from closest mash hits')
        d_qry_to_refs = dict()
        for gid, (rep_gids, contigs_removed) in tqdm(d_gid_to_jobs.items()):
            d_qry_to_refs[gid] = rep_gids
        ref_gids_q = [(x, rep_gid_dir) for x in sorted(ref_gids)]

        log('Caching reference genomes to disk')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            list(tqdm(pool.imap_unordered(genome_cache_worker, ref_gids_q), total=len(ref_gids_q), smoothing=0.05))

        log('Running FastANI')
        for qry_gid, ref_gids in tqdm(d_qry_to_refs.items()):
            path_out = os.path.join(get_gid_root(qry_gid), f'fastani_congruence_c{congruence}_t{target_pct}__k{mash_k}_s{mash_s}.h5')
            run_fastani_on_jobs(qry_gid, ref_gids, tmp_dir, rep_gid_dir, path_out)

    return

def collect_result_worker(job):
    qry_gid, congruence, target_pct, mash_k, mash_s = job
    path_h5 = os.path.join(get_gid_root(qry_gid), f'fastani_congruence_c{congruence}_t{target_pct}__k{mash_k}_s{mash_s}.h5')

    if not os.path.isfile(path_h5):
        raise Exception(f'?? {qry_gid}')

    df = pd.read_hdf(path_h5)

    return df

