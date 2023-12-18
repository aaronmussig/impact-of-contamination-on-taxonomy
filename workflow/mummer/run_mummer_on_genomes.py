import multiprocessing as mp
import os
import tempfile
from typing import FrozenSet

import pandas as pd
from luigi import LocalTarget
from magna.util.disk import move_file, get_file_size_fmt, copy_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, DIR_CACHE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_mummer_on_genome import run_nucmer, COORD_HEADER
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait


class RunMummerOnGenomes(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running Mummer on genomes', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Generating queue')
        ref_gids = df_css.index
        rep_gids = frozenset(df_meta[df_meta['gtdb_representative'] == 't'].index)

        if DIR_CACHE is None:
            raise Exception('needed!')

        cache_dir = os.path.join(DIR_CACHE, 'genomes')

        os.makedirs(cache_dir, exist_ok=True)


        # Run each of the batches
        log('Submitting batches to RQ...')
        queue = [(x, ) for x in ref_gids]
        rq_and_wait(job_id=self.fqn, fn=run_on_gid, q_args=queue, queue_name=self.fqn)

        #
        # for r_gid in tqdm(ref_gids):
        #     # if r_gid != 'GCA_002170165.1':
        #     #     continue
        #     run_on_gid(r_gid)

        # if not DEBUG:
        #     self.write_sentinel()
        return


def cache_gid_worker(job):
    gid, cache_root = job
    from_path = os.path.join(get_gid_root(gid), f'{gid}.fna')
    to_path = os.path.join(get_gid_root(gid, cache_root), f'{gid}.fna')
    os.makedirs(os.path.dirname(to_path), exist_ok=True)

    if not os.path.isfile(to_path):
        copy_file(from_path, to_path, checksum=True)
    return


def run_nucmer_worker(job):
    q_gid, r_gid, path_tmp, cache_dir = job
    return run_nucmer(q_gid, r_gid, path_tmp, cache_dir)


def run_on_gid(ref_gid: str):
    df_meta = GtdbMetadataR207().output().read()
    qry_gids = frozenset(df_meta[df_meta['gtdb_representative'] == 't'].index)
    del df_meta

    cache_dir = os.path.join(DIR_CACHE, 'genomes')

    if not os.path.isdir(cache_dir):
        log(f'Caching files to temporary disk: {cache_dir}')
        to_cache = tuple([(x, cache_dir) for x in sorted(qry_gids)])
        with mp.Pool(processes=mp.cpu_count()) as pool:
            list(tqdm(pool.imap_unordered(cache_gid_worker, to_cache), total=len(to_cache), smoothing=0.01))

    log(f'Processing: {ref_gid}')
    path_out = os.path.join(get_gid_root(ref_gid), 'nucmer_to_rep_genomes.h5')
    if os.path.isfile(path_out):
        return

    # cache the ref gid
    cache_gid_worker((ref_gid, cache_dir))

    qry_gids = tuple(sorted(qry_gids - {ref_gid}))

    with tempfile.TemporaryDirectory() as path_tmp:
        queue = [(q_gid, ref_gid, path_tmp, cache_dir) for q_gid in qry_gids]

        if DEBUG:
            queue = queue[:50]

        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(run_nucmer_worker, queue), total=len(queue), smoothing=0.01))

    rows = list()
    [rows.extend(x) for x in results]
    df = pd.DataFrame(rows, columns=COORD_HEADER)

    if not DEBUG:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = os.path.join(tmp_dir, 'tmp.h5')
            df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')
            log(f'Copying {get_file_size_fmt(tmp_path)} to {path_out}')
            move_file(tmp_path, path_out, checksum=True)
    return
