import multiprocessing as mp
import os
import random
import tempfile
from typing import Dict, Set

import pandas as pd
from Bio import SeqIO
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.bootstrap.select_genomes_for_bootstrapping import SelectGenomesForBootstrapping
from workflow.config import DEBUG, PCT_VALUES, DIR_OUT_PPLACER_RANDOM
from workflow.config import DIR_OUT_BATCH
from workflow.config import DIR_OUT_BOOTSTRAP, DIR_OUT_SENTINEL
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.method.randomly_select_contigs_up_to_pct import randomly_select_contigs_up_to_pct
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.pplacer_random.select_gids_for_pplacer import SelectGidsForPplacerRandom
from workflow.util.collection import iter_batches
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait, rq_wait_for_queue_empty

N_CPUS = mp.cpu_count() // 2

PCT = 50


def worker(job):
    rep, gid_i, gid = job

    gid_root = get_gid_root(gid)
    path_fna = os.path.join(gid_root, f'{gid}.fna')

    d_fna = read_fasta(path_fna)

    for pct in reversed(PCT_VALUES):
        contigs_to_keep, contigs_to_remove = randomly_select_contigs_up_to_pct(d_fna, pct)
        if len(contigs_to_keep) > 0:
            return rep, gid_i, gid, pct, contigs_to_remove

    raise Exception(f'No contigs to keep for {rep} {gid_i} {gid}')




class SelectContigsForPplacer(LuigiTask):

    def requires(self):
        return {
            'bootstrap_gids': SelectGidsForPplacerRandom(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_PPLACER_RANDOM, 'randomly_selected_contigs.h5'))

    def run(self):
        log('Randomly selecting contigs for all genomes at various repeats/samplings', title=True)
        self.make_output_dirs()

        log('Loading bootstrap genomes')
        df_bootstrap = self.input()['bootstrap_gids'].read() if not DEBUG else self.input()['bootstrap_gids'].read_cached()

        log('Creating queue')
        queue = list()
        for rep, row in df_bootstrap.iterrows():
            genomes = row['genomes'].split('|')
            for gid_i, genome in enumerate(genomes):
                queue.append((rep, gid_i, genome))

        if DEBUG:
            queue = queue[:3]

        log(f'Processing queue ({len(queue):,})')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log('Saving to disk')
        rows = list()
        for rep, gid_i, gid, pct, contigs_to_remove in results:
            rows.append({
                'rep': rep,
                'gid_i': gid_i,
                'gid': gid,
                'pct': pct,
                'contigs_to_remove': '|'.join(sorted(contigs_to_remove)),
            })

        log('Creating dataframe')
        df = pd.DataFrame(rows)

        if not DEBUG:
            log('Saving to disk')
            self.save_hdf(df)


