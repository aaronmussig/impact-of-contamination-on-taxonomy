import multiprocessing as mp
import os
import subprocess
import tempfile

from Bio import SeqIO
from luigi import LocalTarget
from magna.util.disk import copy_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, \
    DIR_OUT_BLAST_REP_DB
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def collect_fna_worker(gid):
    fna_path = os.path.join(get_gid_root(gid), f'{gid}.fna')
    out = list()
    for record in SeqIO.parse(fna_path, 'fasta'):
        out.append((gid, record.id, str(record.seq)))
    return out


class MakeRepBlastDb(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def load_reps(self):
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        return sorted(set(df_meta[df_meta['gtdb_representative'] == 't'].index))

    def run(self):
        log('Creating blast database from GTDB R207 representative genomes', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        rep_gids = self.load_reps()

        batches = list(iter_batches(rep_gids, 1000))
        gids_seen = set()
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_fna_tmp = os.path.join(tmp_dir, 'concat.fna')
            log(f'Writing to: {path_fna_tmp}')
            with open(path_fna_tmp, 'w') as f:
                for batch in tqdm(batches):
                    with mp.Pool(processes=min(mp.cpu_count(), 40)) as pool:
                        for results in pool.imap_unordered(collect_fna_worker, batch):
                            for gid, contig, seq in results:
                                f.write(f'>{gid}_{contig}\n{seq}\n')
                                gids_seen.add(gid)

            log(f'Processed: {len(gids_seen):,}')
            print(f'Concat FNA: {get_file_size_fmt(path_fna_tmp)}')

            db_out_dir = os.path.join(tmp_dir, 'blast_db')
            os.makedirs(db_out_dir)
            path_blast_db_tmp = os.path.join(db_out_dir, 'blastdb')

            cmd = [
                'makeblastdb',
                '-dbtype', 'nucl',
                '-in', path_fna_tmp,
                '-title', 'GTDB_R207_REPS',
                '-parse_seqids',
                '-out', path_blast_db_tmp,
                '-max_file_sz', '4GB',
                '-hash_index'
            ]

            subprocess.check_output(cmd)

            log(f'Copying to: {DIR_OUT_BLAST_REP_DB}')
            for file in os.listdir(db_out_dir):
                copy_file(os.path.join(db_out_dir, file), os.path.join(DIR_OUT_BLAST_REP_DB, file), checksum=True)

        return
