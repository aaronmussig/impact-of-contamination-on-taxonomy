import multiprocessing as mp
import os
import shutil
import tempfile
from collections import defaultdict
from typing import Dict, Set, List

import pandas as pd
from chiton.fastani import fastani
from magna.util.disk import get_file_size_fmt, move_file
from tqdm import tqdm

from workflow.config import DIR_OUT_V2_FASTANI_REP_TO_CLOSEST_REP, DIR_OUT_BATCH, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask, LocalTargetTsv, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait, rq_wait_for_queue_empty


class V2FastAniRepToClosestRepCreateJobs(LuigiTask):

    @property
    def root_dir(self):
        return DIR_OUT_V2_FASTANI_REP_TO_CLOSEST_REP

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, 'create_jobs_results.h5'))

    def run(self):
        log(f'V2FastAniRepToClosestRepCreateJobs()', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        df_reps = df_meta[(df_meta['gtdb_representative'] == 't')]
        rep_gids = set(df_reps.index)
        log(f'Found {len(rep_gids):,} rep gids')

        log('Creating rank to genome id mapping (reps only)')
        d_rank_to_gids = get_rank_to_genome_id(df_reps)

        # Create the comparison list
        d_ref_to_query = dict()
        for row in df_reps.itertuples():
            for rank in reversed(row.gtdb_taxonomy.split(';')):
                n_reps_in_rank = len(d_rank_to_gids[rank])
                if n_reps_in_rank > 1:
                    d_ref_to_query[row.Index] = d_rank_to_gids[rank]
                    break

        n_comparisons = sum(len(x) for x in d_ref_to_query.values())
        log(f'Creating jobs for {n_comparisons:,} comparisons')

        batch_root = os.path.join(DIR_OUT_BATCH, f'{self.fqn}')
        job_id = self.fqn
        if os.path.isdir(batch_root):
            batch_paths = list()
            for ref, qry in tqdm(d_ref_to_query.items()):
                path = os.path.join(batch_root, f'{ref}.tsv')
                batch_paths.append((path,))
        else:
            os.makedirs(batch_root, exist_ok=True)
            batch_paths = list()
            for ref, qry in tqdm(d_ref_to_query.items()):
                path = os.path.join(batch_root, f'{ref}.tsv')
                with open(path, 'w') as f:
                    for gid in sorted(qry):
                        f.write(f'{gid}\n')
                batch_paths.append((path,))

            log('Submitting batches to RQ...')
            log(f'Job ID: {job_id}')
            rq_and_wait(job_id=job_id, fn=run_on_batch_file, q_args=batch_paths, queue_name=job_id)

        log('Waiting until RQ queue is empty...')
        rq_wait_for_queue_empty(job_id)

        # Collect the results and clean-up the batch directory
        log('Collecting results...')
        dfs: List[pd.DataFrame] = list()
        path_queue = list()
        for path in tqdm(batch_paths):
            path = f'{path[0][0:-4]}.h5'
            path_queue.append(path)
        with mp.Pool(processes=mp.cpu_count()) as pool:
            for df in tqdm(pool.imap(pd.read_hdf, path_queue)):
                dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        log('Writing results...')
        self.save_hdf(df)

        log('Cleaning up batch directory...')
        print(batch_root)
        # shutil.rmtree(batch_root)

        log('Done.')
        return


def get_rank_to_genome_id(df_reps: pd.DataFrame) -> Dict[str, Set[str]]:
    out = defaultdict(set)
    for row in df_reps.itertuples():
        gid = row.Index
        for rank in row.gtdb_taxonomy.split(';'):
            out[rank].add(gid)
    return dict(out)


def run_on_batch_file(path):
    ref_gids = {os.path.basename(path).replace('.tsv', '')}
    assert (len(ref_gids) == 1)
    query_gids = set()
    with open(path) as f:
        for line in f.readlines():
            query_gids.add(line.strip())

    df = run_fastani_on_jobs(query_gids, ref_gids)
    out_path = f'{path[0:-4]}.h5'
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, 'out.h5')
        df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')

        log(f'Copying {get_file_size_fmt(tmp_path)} to {out_path}')
        move_file(tmp_path, out_path, checksum=True)
    return


def run_fastani_on_jobs(ref_gids: Set[str], qry_gids: Set[str]):
    ref_gid_paths = {x: os.path.join(get_gid_root(x), f'{x}.fna') for x in ref_gids}
    qry_gid_paths = {x: os.path.join(get_gid_root(x), f'{x}.fna') for x in qry_gids}

    out = list()

    ani = fastani(query=list(qry_gid_paths.values()),
                  reference=list(ref_gid_paths.values()),
                  cpus=mp.cpu_count(),
                  single_execution=False,
                  bidirectional=True,
                  show_progress=True)

    # Prepare the ANI results for output
    d_ani = ani.as_dict()
    for qry_gid in qry_gid_paths.keys():
        q_key = qry_gid_paths[qry_gid]
        for ref_gid in ref_gid_paths.keys():
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
            })

    df = pd.DataFrame(out)
    return df
