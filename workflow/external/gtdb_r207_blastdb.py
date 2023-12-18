import multiprocessing as mp
import os
import subprocess
import tempfile

from luigi import LocalTarget
from magna.util.disk import copy_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DIR_OUT_EXTERNAL_R207_BLAST_DB, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class GtdbR207BlastDb(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_EXTERNAL_R207_BLAST_DB, 'gtdb_r207_blast.db.ndb'))

    def run(self):
        log('Creating GTDB R207 BLAST DB', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Creating queue')
        queue = sorted(df_meta.index)
        del df_meta

        out_dir = os.path.dirname(self.output().path)

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_fna_local = os.path.join(tmp_dir, 'gtdb_r207.fna')
            path_db_out = os.path.join(tmp_dir, 'gtdb_r207_blast.db')
            log(f'Creating temporary FNA: {path_fna_local}')

            # with mp.Pool(processes=min(40, mp.cpu_count())) as pool:
            with open(path_fna_local, 'w') as f:
                for gid in tqdm(queue):
                    _, seq_data = collect_fna_worker(gid)
                    for contig, seq in seq_data:
                        f.write(f'>{gid}__{contig}\n{seq}\n')

            log('Creating blast database')
            cmd = [
                'makeblastdb',
                '-in', path_fna_local,
                '-dbtype', 'nucl',
                '-blastdb_version', '5',
                '-out', path_db_out,
            ]
            subprocess.run(cmd, check=True)

            log('Copying to output')
            for file_name in os.listdir(tmp_dir):
                file_path = os.path.join(tmp_dir, file_name)
                if file_path != path_fna_local:
                    print(f'Copying {file_path} to {out_dir} ({get_file_size_fmt(file_path)})')
                    if not DEBUG:
                        copy_file(file_path, os.path.join(out_dir, file_name), checksum=True)

        return


def collect_fna_worker(gid):
    genome = Genome(gid)
    return gid, list(genome.d_fna.items())
