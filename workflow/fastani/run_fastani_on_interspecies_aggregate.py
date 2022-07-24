import multiprocessing as mp
import os
import tempfile
from collections import defaultdict

import pandas as pd
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTANI
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani.run_fastani_on_interspecies_ani import RunFastAniOnInterspeciesAni
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root

N_CPUS = mp.cpu_count()


class RunFastAniOnInterspeciesAniAggregate(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            '_ani': RunFastAniOnInterspeciesAni()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_FASTANI, 'fastani_interspecies_agg.h5'))

    def run(self):
        log('Running FastANI interspecies ANI for each genus (aggregate)', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        df_meta = df_meta[df_meta['gtdb_representative'] == 't']

        log('Mapping genus to species representatives')
        d_genus_to_species_reps = defaultdict(set)
        for gid, row in tqdm(df_meta.iterrows(), total=len(df_meta)):
            d_genus_to_species_reps[row['genus']].add(gid)

        log('Creating queue')
        genera = [k for k, v in d_genus_to_species_reps.items() if len(v) > 1]

        if DEBUG:
            genera = genera[0:10]

        log('Processing queue')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(worker, genera), total=len(genera)))

        log('Writing results')
        df = pd.concat(results)
        self.save_hdf(df, data_columns=['genus'])


def worker(genus):
    df = pd.read_hdf(os.path.join(DIR_OUT_FASTANI, 'interspecies_ani', f'{genus}.h5'))
    df['genus'] = genus
    return df[['genus', 'query', 'reference', 'ani', 'af']]


def copy_genome_worker(job):
    gid, tmp_dir = job
    gid_path_srv = os.path.join(get_gid_root(gid), f'{gid}.fna')
    gid_path_tmp = os.path.join(tmp_dir, f'{gid}.fna')
    if not os.path.isfile(gid_path_tmp):
        copy_file(gid_path_srv, gid_path_tmp)


def fastani_worker(genus, sp_reps, tmp_dir):
    # Read the FASTA file
    path_results_out = os.path.join(DIR_OUT_FASTANI, 'interspecies_ani', f'{genus}.h5')

    # Stop early if this already exists
    if os.path.isfile(path_results_out):
        return

    # Convert the rep set into their paths
    d_gid_to_path = dict()
    d_path_to_gid = dict()
    gids = list()
    for gid in sp_reps:
        d_gid_to_path[gid] = os.path.join(tmp_dir, f'{gid}.fna')
        d_path_to_gid[d_gid_to_path[gid]] = gid
        gids.append(gid)

    # Prepare for FastANI
    ani = fastani(query=list(d_gid_to_path.values()),
                  reference=list(d_gid_to_path.values()),
                  cpus=N_CPUS,
                  single_execution=False,
                  bidirectional=True,
                  show_progress=False)

    # Prepare the ANI results for output
    cur_results = list()
    d_ani = ani.as_dict()

    for i in range(len(gids)):
        for j in range(i + 1):
            qry_gid, ref_gid = gids[i], gids[j]
            qry_path, ref_path = d_gid_to_path[qry_gid], d_gid_to_path[ref_gid]

            qvr = d_ani[qry_path][ref_path]
            rvq = d_ani[ref_path][qry_path]
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
    df.to_hdf(path_results_tmp, key='root', format='table', complevel=9, complib='blosc:lz4hc')
    move_file(path_results_tmp, path_results_out, checksum=True)
    return


def batch_worker(batch_path):
    log(f'Reading batch file: {batch_path}')
    queue = list()
    gids_to_cache = set()
    with open(batch_path) as f:
        for line in f.readlines():
            genus, sp_reps = line.strip().split('\t')
            sp_reps = frozenset(sp_reps.split('|'))
            gids_to_cache.update(sp_reps)
            queue.append((genus, sp_reps))
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
        for genus, sp_reps in tqdm(queue):
            fastani_worker(genus, sp_reps, tmp_dir)
    return
