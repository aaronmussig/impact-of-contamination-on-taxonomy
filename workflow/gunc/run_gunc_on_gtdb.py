import multiprocessing as mp
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Tuple, List

from luigi import LocalTarget
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import PATH_GUNC_GTDB_REF_DB, DIR_OUT_SENTINEL, DIR_OUT_BATCH, DEBUG, DIR_CACHE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.genome.symlink import GenomeSymlink
from workflow.model.luigi import LuigiTask
from workflow.util.collection import iter_batches
from workflow.util.gunc import move_gunc_directory, gunc_are_all_files_present
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait


def batch_worker(batch_path):
    # Cache the database to disk if required
    if DIR_CACHE is not None:
        path_gunc_db = os.path.join(DIR_CACHE, 'gunc_db_gtdb95.dmnd')
        os.makedirs(DIR_CACHE, exist_ok=True)
        if not os.path.isfile(path_gunc_db):
            log(f'Copying {PATH_GUNC_GTDB_REF_DB} to {path_gunc_db}')
            copy_file(PATH_GUNC_GTDB_REF_DB, path_gunc_db, checksum=True)
            log('Done.')
    else:
        path_gunc_db = PATH_GUNC_GTDB_REF_DB

    # Get the gid mapping for each batch
    d_gid_mapping = dict()
    with open(batch_path) as f:
        for line in f.readlines():
            path = line.strip()
            prodigal_dir = Path(os.path.dirname(path))
            gunc_dir = os.path.join(prodigal_dir.parent.absolute(), 'gunc_r95')
            gid = os.path.basename(path).replace('.faa', '')
            d_gid_mapping[gid] = gunc_dir

    # Run GUNC into a temporary directory
    with tempfile.TemporaryDirectory() as f_out, tempfile.TemporaryDirectory() as f_tmp:
        cmd = ['gunc', 'run', '--db_file', path_gunc_db,
               '--file_suffix', '.faa', '--gene_calls',
               '--input_file', batch_path, '--threads',
               str(mp.cpu_count()), '--detailed_output',
               '--contig_taxonomy_output', '--out_dir', f_out,
               '--temp_dir', f_tmp]
        subprocess.check_output(cmd)

        # Move the output files to the output directory
        move_gunc_directory(f_out, d_gid_mapping)
        return


def create_batches(faa_paths, root_dir, size) -> List[Tuple[str]]:
    """Partition the FAA files into batches."""
    log(f'Creating batches of size={size:,} in {root_dir} for {len(faa_paths):,} files')
    out = list()
    for i, cur_batch in enumerate(iter_batches(sorted(faa_paths), size)):
        cur_path = os.path.join(root_dir, f'batch_{i}.txt')
        with open(cur_path, 'w') as f:
            f.write('\n'.join(cur_batch) + '\n')
        out.append((cur_path,))
    return out


def get_path_worker(gid):
    gid_root = get_gid_root(gid)
    gunc_dir_r95 = os.path.join(gid_root, 'gunc_r95')
    if not gunc_are_all_files_present(gid, gunc_dir_r95, db='gtdb_95'):
        return os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    return None


class RunGuncOnGtdbR95UsingR207(LuigiTask):
    """Runs GUNC on the GTDB R207 genomes using RQ."""

    BATCH_SIZE = 1000 if not DEBUG else 2

    def requires(self):
        return {
            'genomes': GenomeSymlink(),
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running GUNC on GTDB R95 using R207...', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        gids = sorted(set(df_meta.index))

        # Get the queue of gid, FAA to process
        log('Getting called gene paths')
        if DEBUG:
            faa_paths = [get_path_worker(gid) for gid in gids]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                faa_paths = list(tqdm(pool.imap_unordered(get_path_worker, gids), total=len(gids)))
        faa_paths = [x for x in faa_paths if x is not None]
        log(f'Found {len(faa_paths):,} FAA files to process.')

        log('Creating batches...')
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        os.makedirs(batch_dir, exist_ok=True)
        batch_paths = create_batches(faa_paths, root_dir=batch_dir, size=self.BATCH_SIZE)
        log(f'Created {len(batch_paths):,} batches ~{self.BATCH_SIZE:,} genomes each.')

        # Run each of the batches
        log('Submitting batches to RQ...')
        rq_and_wait(job_id=self.fqn, fn=batch_worker, q_args=batch_paths, queue_name=self.fqn)

        if not DEBUG:
            self.write_sentinel()
