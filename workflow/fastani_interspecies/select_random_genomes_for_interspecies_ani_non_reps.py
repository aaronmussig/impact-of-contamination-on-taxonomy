import multiprocessing as mp
import os
import tempfile
from collections import defaultdict

import pandas as pd
from chiton.fastani import fastani
from luigi import LocalTarget
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, DIR_OUT_FASTANI_INTER
from workflow.config import DIR_OUT_BATCH
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani.remove_gunc_failed_contigs_by_contamination_sp_cluster import \
    RemoveGuncFailedContigsByContaminationSpCluster
from workflow.model.luigi import LuigiTask, LocalTargetTsv, LocalTargetHdf5
from workflow.util.collection import invert_dict
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty
import random

N_CPUS = mp.cpu_count()

N_RANDOM = 10


class SelectRandomGenomesForInterspeciesAniNonReps(LuigiTask):
    """
    Run FastANI on all sp cluster failed genomes to all species reps within the genus
    """

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_FASTANI_INTER, 'random_genomes_for_interspecies.h5'))

    def run(self):
        log('Selecting random genomes for species ani non reps', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        df_meta_non_reps = df_meta[df_meta['gtdb_representative'] == 'f']

        d_genus_to_non_reps = defaultdict(set)
        for gid, row in tqdm(df_meta_non_reps.iterrows(), total=len(df_meta_non_reps)):
            d_genus_to_non_reps[row['genus']].add(gid)

        rows = list()
        for genus, gids in d_genus_to_non_reps.items():
            if N_RANDOM >= len(gids):
                gids_selected = gids
            else:
                gids_selected = set(random.sample(gids, k=N_RANDOM))

            rows.append({
                'genus': genus,
                'gids': '|'.join(sorted(gids_selected))
            })

        df = pd.DataFrame(rows)

        self.save_hdf(df, index='genus')
        return



