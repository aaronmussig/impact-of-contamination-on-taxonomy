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
    prodigal_dir = os.path.join(root, 'prodigal')
    tigr_th = os.path.join(prodigal_dir, f'{gid}_tigrfam_tophit.tsv')
    pfam_th = os.path.join(prodigal_dir, f'{gid}_pfam_tophit.tsv')
    os.makedirs(prodigal_dir, exist_ok=True)

    # Get the remote targets
    r207_root = get_gid_r207_root(gid)

    if not os.path.exists(tigr_th):
        old_tigr_th = os.path.join(r207_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv')
        if not os.path.exists(old_tigr_th):
            raise Exception(f'Could not find {old_tigr_th}')
        os.symlink(old_tigr_th, tigr_th)

    if not os.path.exists(pfam_th):
        old_pfam_th = os.path.join(r207_root, 'prodigal', f'{gid}_pfam_tophit.tsv')
        if not os.path.exists(old_pfam_th):
            raise Exception(f'Could not find {old_pfam_th}')
        os.symlink(old_pfam_th, pfam_th)

    return


class GenomeSymlinkMarker(LuigiTask):

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

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        queue = sorted(set(df_meta.index))
        log(f'Processing {len(queue):,} genomes')

        if DEBUG:
            [worker(x) for x in queue]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))
            self.write_sentinel()
