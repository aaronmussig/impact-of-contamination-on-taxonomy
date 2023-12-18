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

from workflow.config import DEBUG, DIR_OUT_BATCH, DIR_CACHE, DIR_OUT_FASTANI_CONTIG_SPLIT_RANDOM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani_contig_split.random_a_run_mash import FastAniContigSplitRandomRunMash
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_wait_for_queue_empty, rq_and_wait


class FastAniContigSplitRandomRunFastAni(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)
    replicate_id = luigi.IntParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTANI_CONTIG_SPLIT_RANDOM, f'r_{self.replicate_id}')

    def requires(self):
        return {
            'mash': FastAniContigSplitRandomRunMash(
                target_pct=self.target_pct,
                replicate_id=self.replicate_id,
                mash_k=self.mash_k,
                mash_s=self.mash_s
            ),
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'closest_reps': Closest100GenomesToRepresentative()
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                self.root_dir,
                f'ani_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}_r{self.replicate_id}.h5'
            ))

    def run(self):
        log(f'FastAniContigSplitRandomRunFastAni(p={self.target_pct}, r={self.replicate_id})', title=True)
        self.make_output_dirs()

        batch_size = 25

        batch_root = os.path.join(DIR_OUT_BATCH,
                                  f'{self.fqn}_k{self.mash_k}_s{self.mash_s}___p{self.target_pct}_r{self.replicate_id}')
        # batch_root = '/tmp/batch' if DEBUG else batch_root

        if not os.path.isdir(batch_root):
            log('Creating batch paths')
            os.makedirs(batch_root, exist_ok=True)
            batch_paths = self.create_batch_files(batch_root, batch_size)
            rq_submit = True
        else:
            log('Loading batch paths')
            batch_paths = self.load_batch_files(batch_root)
            rq_submit = False

        if DEBUG:
            for batch_path in batch_paths:
                run_on_batch_file(batch_path, self.replicate_id, self.target_pct, self.mash_k, self.mash_s)

        if rq_submit:
            queue = [(x, self.replicate_id, self.target_pct, self.mash_k, self.mash_s) for x in batch_paths]
            log('Submitting batches to RQ...')
            job_id = f'{self.fqn}_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}_r{self.replicate_id}'
            log(f'Job ID: {job_id}')
            rq_and_wait(job_id=job_id, fn=run_on_batch_file, q_args=queue, queue_name=job_id)
            # self.write_sentinel()

            # Run each of the batches
            log('Waiting until RQ queue is empty...')
            rq_wait_for_queue_empty(job_id)

        # Collecting the results
        log('Collecting results')
        log('Reading batch files')
        queue = list()
        for job, _ in (read_batch_file(x) for x in tqdm(batch_paths)):
            for qry_gid in job.keys():
                queue.append((qry_gid, self.target_pct, self.replicate_id, self.mash_k, self.mash_s))

            if DEBUG and len(queue) > 5:
                break

        log('Reading dataframes')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            dfs = list(tqdm(pool.imap_unordered(collect_result_worker, queue), total=len(queue)))

        log('Concat dfs')
        df = pd.concat(dfs, ignore_index=True)
        df.sort_values(by=['query', 'ani', 'af', 'ref'], ascending=[True, False, False, True], inplace=True,
                       ignore_index=True)
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

        log('Reading MAXCSS file')
        df_css = self.input()['max_css'].maybe_read_cached()

        d_gid_to_seed = {k: i + self.replicate_id for i, k in enumerate(df_css.index)}

        log('Allocating genomes into batches')
        d_batch_id_to_gids = defaultdict(set)
        for batch_id, batch_gids in enumerate(iter_batches(sorted(df_css.index), n=batch_size)):
            for batch_gid in batch_gids:
                d_batch_id_to_gids[batch_id].add(batch_gid)

        d_qry_to_refs = defaultdict(set)
        log('Adding genomes from the expected taxa')
        for qry_gid in tqdm(df_css.index):
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
        for qry_gid in tqdm(df_css.index):
            cur_sp_rep_gid = df_meta.loc[qry_gid, 'gtdb_genome_representative'][3:]
            closest_reps = frozenset(df_closest_reps.loc[cur_sp_rep_gid, 'closest_representatives'].split('|'))
            d_qry_to_refs[qry_gid].update(closest_reps)

        log('Loading Mash scores')
        df_mash = self.input()['mash'].maybe_read_cached()
        for row in tqdm(df_mash.itertuples(), total=len(df_mash)):
            # Note; This :15 makes sure that the suffix _C is removed.
            # We will want to run each half against both sets of genomes.
            d_qry_to_refs[row.query[:15]].add(row.ref)

        log('Creating batch files')
        os.makedirs(batch_root, exist_ok=True)
        batch_paths = list()
        for batch_id, batch_gids in d_batch_id_to_gids.items():
            batch_path = os.path.join(batch_root, f'{str(batch_id)}.tsv')
            batch_paths.append(batch_path)
            with open(batch_path, 'w') as f:
                for batch_gid in batch_gids:
                    gid_seed = d_gid_to_seed[batch_gid]
                    rep_gids_str = '|'.join(d_qry_to_refs[batch_gid])
                    if len(rep_gids_str) < 2:
                        raise Exception(f'{batch_gid} has no mash hits')
                    else:
                        f.write(f'{batch_gid}\t{gid_seed}\t{rep_gids_str}\n')
        return batch_paths

    def load_batch_files(self, batch_root):
        out = list()
        for file in os.listdir(batch_root):
            if file.endswith('.tsv'):
                out.append(os.path.join(batch_root, file))
        return out


def run_fastani_on_jobs(ref_gids: Set[str], tmp_ref_dir, path_out: str, keep_fna, disc_fna):
    if os.path.exists(path_out):
        return

    ref_gid_paths = {x: os.path.join(get_gid_root(x, root=tmp_ref_dir), f'{x}.fna') for x in ref_gids}

    out = list()
    for path_fna_tmp in (keep_fna, disc_fna):

        ani = fastani(query=path_fna_tmp,
                      reference=list(ref_gid_paths.values()),
                      cpus=mp.cpu_count(),
                      single_execution=False,
                      bidirectional=True,
                      show_progress=False)

        # Prepare the ANI results for output
        d_ani = ani.as_dict()
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
                'query': os.path.basename(path_fna_tmp).replace('.fna', ''),
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
    gid, pct_to_remove, seed, tmp_dir = job

    genome = Genome(gid)
    d_contig_to_seq_keep, d_contig_to_seq_disc = genome.split_fna_by_pct_random(pct=pct_to_remove, seed=seed)

    path_fna_keep = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}.fna')
    os.makedirs(os.path.dirname(path_fna_keep), exist_ok=True)

    path_fna_disc = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}_C.fna')
    os.makedirs(os.path.dirname(path_fna_disc), exist_ok=True)

    # Write the modified FNA files to disk
    with open(path_fna_keep, 'w') as f:
        for contig, seq in d_contig_to_seq_keep.items():
            f.write(f'>{contig}\n{seq}\n')

    with open(path_fna_disc, 'w') as f:
        for contig, seq in d_contig_to_seq_disc.items():
            f.write(f'>{contig}\n{seq}\n')

    return gid, path_fna_keep, path_fna_disc


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
            gid, seed, rep_gids_str = line.strip().split('\t')
            rep_gids = frozenset(rep_gids_str.split('|'))
            d_gid_to_jobs[gid] = (int(seed), rep_gids)
            ref_gids.update(rep_gids)
    return d_gid_to_jobs, ref_gids


def run_on_batch_file(batch_path, replicate_id, target_pct, mash_k, mash_s):
    log(f'Processing batch file: {batch_path}')
    d_gid_to_jobs, ref_gids = read_batch_file(batch_path)

    rep_gid_dir = os.path.join(DIR_CACHE, 'genomes')

    with tempfile.TemporaryDirectory() as tmp_dir:

        log('Getting the pruned genome file for all failed')
        queue = list()
        queued_gids = set()
        for gid, (seed, rep_gids) in d_gid_to_jobs.items():
            queue.append((
                gid,
                target_pct,
                seed,
                tmp_dir,
            ))
            queued_gids.add(gid)

        log('Processing queue')
        if DEBUG:
            genome_paths = [write_genome_fna_worker(x) for x in tqdm(queue)]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                genome_paths = list(
                    tqdm(pool.imap_unordered(write_genome_fna_worker, queue), total=len(queue), smoothing=0.05))
        d_gid_to_paths = {x[0]: (x[1], x[2]) for x in genome_paths}

        log('Generating queue from closest mash hits')
        d_qry_to_refs = dict()
        for gid, (seed, rep_gids) in tqdm(d_gid_to_jobs.items()):
            keep_fna, disc_fna = d_gid_to_paths[gid]
            d_qry_to_refs[gid] = (rep_gids, keep_fna, disc_fna)
        ref_gids_q = [(x, rep_gid_dir) for x in sorted(ref_gids)]

        log('Caching reference genomes to disk')
        if DEBUG:
            print('skipping')
            # [genome_cache_worker(x) for x in tqdm(ref_gids_q)]
        else:
            with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
                list(tqdm(pool.imap_unordered(genome_cache_worker, ref_gids_q), total=len(ref_gids_q), smoothing=0.05))

        log('Running FastANI')
        for qry_gid, (rep_gids, keep_fna, disc_fna) in tqdm(d_qry_to_refs.items()):
            path_out = os.path.join(get_gid_root(qry_gid),
                                    f'fastani_contig_split_random_t{target_pct}__k{mash_k}_s{mash_s}_r{replicate_id}.h5')
            run_fastani_on_jobs(rep_gids, rep_gid_dir, path_out, keep_fna, disc_fna)

    return


def collect_result_worker(job):
    qry_gid, target_pct, replicate_id, mash_k, mash_s = job
    path_h5 = os.path.join(get_gid_root(qry_gid),
                           f'fastani_contig_split_random_t{target_pct}__k{mash_k}_s{mash_s}_r{replicate_id}.h5')

    if not os.path.isfile(path_h5):
        raise Exception(f'?? {qry_gid}')

    df = pd.read_hdf(path_h5)

    return df[df['ani'] > 0]
