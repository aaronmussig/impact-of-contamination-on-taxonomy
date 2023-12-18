import multiprocessing as mp
import os
import subprocess
import tempfile

import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_BLAST
from workflow.external.gtdb_r207_blastdb import GtdbR207BlastDb
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class BlastGenomeAgainstR207Db(LuigiTask):

    def requires(self):
        return {
            'blast': GtdbR207BlastDb(),
            'css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_BLAST, f'blast_genome_against_207_db_tmp.h5'))

    def run(self):
        log(f'BlastGenomeAgainstR207Db', title=True)

        log('Loading MaxCSS file')
        df_css = self.input()['css'].maybe_read_cached()

        log('Creating query file')
        queue = sorted(set(df_css.index))
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(read_genome_fna, queue), total=len(queue)))
        #
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_query = os.path.join(tmp_dir, 'query.fna')
            log(f'Writing query file to {path_query}')
            with open(path_query, 'w') as f:
                for gid, result in results:
                    if gid != 'GCA_900759445.1':
                        continue

                    for contig, seq in result.items():
                        f.write(f'>{gid}__{contig}\n{seq}\n')

            # Free memory
            del results
            del df_css
            del queue

            path_blast_tmp = '/tmp/am/blast_results.out'
            os.makedirs('/tmp/am', exist_ok=True)

            log('Running blast')
            cmd = [
                'blastn',
                '-query',
                path_query,
                '-db',
                self.input()['blast'].path[:-4],
                '-num_threads',
                str(mp.cpu_count()),
                '-outfmt',
                '7',
                '-out',
                path_blast_tmp,
            ]
            log(' '.join(cmd))

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()

            if proc.returncode != 0:
                raise Exception(f'Error running blast: {stderr}')

            log('Parsing results')
            df = parse_blast_output(path_blast_tmp)

        log('Saving results')
        self.save_hdf(df)

        return


def read_genome_fna(gid):
    genome = Genome(gid)
    d_fna = genome.d_fna
    return gid, d_fna


def parse_blast_output(path):
    rows = list()

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')

            q_gid, q_contig = cols[0].split('__')
            s_gid, s_contig = cols[1].split('__')

            rows.append({
                'q_gid': q_gid,
                'q_contig': q_contig,
                'r_gid': s_gid,
                'r_contig': s_contig,
                'pct_identity': float(cols[2]),
                'aln_len': int(cols[3]),
                'mismatches': int(cols[4]),
                'gap_openings': int(cols[5]),
                'q_start': int(cols[6]),
                'q_end': int(cols[7]),
                'r_start': int(cols[8]),
                'r_end': int(cols[9]),
                'e_value': float(cols[10]),
                'bit_score': float(cols[11]),
            })

    df = pd.DataFrame(rows)

    return df
