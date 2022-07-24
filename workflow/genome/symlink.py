import multiprocessing as mp
import os

from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DIR_OUT_SENTINEL, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from workflow.util.paths import get_gid_root, get_gid_r207_root


def worker(gid):
    # Get the destination directory
    root = get_gid_root(gid)

    # Define targets
    fna = os.path.join(root, f'{gid}.fna')
    prodigal_dir = os.path.join(root, 'prodigal')
    prodigal_fna = os.path.join(prodigal_dir, f'{gid}.fna')
    prodigal_faa = os.path.join(prodigal_dir, f'{gid}.faa')
    os.makedirs(prodigal_dir, exist_ok=True)

    # Get the remote targets
    r207_root = get_gid_r207_root(gid)

    if not os.path.exists(fna):
        old_fna = os.path.join(r207_root, f'{os.path.basename(r207_root)}_genomic.fna')
        os.symlink(old_fna, fna)

    if not os.path.exists(prodigal_fna):
        old_prodigal_fna = os.path.join(r207_root, 'prodigal', f'{gid}_protein.fna')
        os.symlink(old_prodigal_fna, prodigal_fna)

    if not os.path.exists(prodigal_faa):
        old_prodigal_fna = os.path.join(r207_root, 'prodigal', f'{gid}_protein.faa')
        os.symlink(old_prodigal_fna, prodigal_faa)


class GenomeSymlink(LuigiTask):
    """Symlink all genome and prodigal files from GTDB R95."""

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def _get_queue(self):
        df_meta = self.input()['meta'].read()
        out = tuple(sorted(set(df_meta.index)))
        log(f'{len(out):,} genomes to process.')
        return out

    def run(self):
        log('Linking all GTDB R207 genomes...', title=True)

        log('Creating output directory')
        self.make_output_dirs()

        log('Loading R207 metadata')
        df = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        queue = sorted(set(df.index))
        log(f'Processing {len(queue):,} genomes')

        if DEBUG:
            [worker(x) for x in queue]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))
            self.write_sentinel()
