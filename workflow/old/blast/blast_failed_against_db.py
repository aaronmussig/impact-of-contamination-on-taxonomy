import multiprocessing as mp
import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from luigi import LocalTarget
from magna.util.disk import get_file_size_fmt, move_file
from tqdm import tqdm

from workflow.old.blast.make_rep_blast_db import MakeRepBlastDb
from workflow.config import DEBUG, DIR_OUT_BLAST
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def gid_fna_worker(gid):
    gid_root = get_gid_root(gid)
    path_fna = os.path.join(gid_root, f'{gid}.fna')

    out = list()
    for record in SeqIO.parse(path_fna, 'fasta'):
        out.append((gid, record.id, str(record.seq)))
    return out


class BlastFailedAgainstDb(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            '_blast': MakeRepBlastDb(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_BLAST, 'fail_blast_against_db.h5'))

    def get_queue(self):
        df_max_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()
        queue = sorted(df_max_css.index)
        log(f'Found: {len(queue):,} genomes to process')
        return queue

    def run(self):
        log('Blasting failed genomes against database', title=True)
        self.make_output_dirs()

        queue = self.get_queue()

        log('Collecting fasta files into one file')
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_concat = os.path.join(tmp_dir, 'concat.fna')
            with open(path_concat, 'w') as f:
                with mp.Pool(processes=mp.cpu_count()) as pool:
                    for results in tqdm(pool.imap_unordered(gid_fna_worker, queue), total=len(queue)):
                        for gid, contig, seq in results:
                            f.write(f'>{gid}_{contig}\n{seq}\n')

            log('Blasting genomes')
            path_blast_tmp = os.path.join(tmp_dir, 'blast.out')
            cmd = [
                'blastn',
                # '-db', os.path.join(DIR_OUT_BLAST_REP_DB, 'blastdb'),
                '-db', '/tmp/am/blastdb',
                '-query', path_concat,
                '-out', path_blast_tmp,
                '-num_threads', str(mp.cpu_count()),
                '-outfmt', '6'
            ]

            subprocess.check_output(cmd)

            header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            df = pd.read_csv(path_blast_tmp, sep='\t', names=header, low_memory=False)

            df_tmp_path = os.path.join(tmp_dir, f'{gid}.h5')
            df.to_hdf(df_tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')

            log(f'Copying {get_file_size_fmt(df_tmp_path)} to {self.output().path}')
            move_file(df_tmp_path, self.output().path, checksum=True)

        return
