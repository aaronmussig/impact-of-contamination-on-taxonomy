import multiprocessing as mp
import os

import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, DIR_OUT_FASTANI_INTER
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani_interspecies.run_fastani_interspecies_to_genus_ani_for_random_non_reps import \
    FastAniInterspeciesToGenusAniForRandomNonReps
from workflow.fastani_interspecies.select_random_genomes_for_interspecies_ani_non_reps import \
    SelectRandomGenomesForInterspeciesAniNonReps
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log

N_CPUS = mp.cpu_count()


class FastAniInterspeciesToGenusAniForRandomNonRepsAgg(LuigiTask):
    """
    Run FastANI on all sp cluster failed genomes to all species reps within the genus
    """

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'selected': SelectRandomGenomesForInterspeciesAniNonReps(),
            '_done': FastAniInterspeciesToGenusAniForRandomNonReps()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running FastANI interspecies genus ANI for random non reps (agg)', title=True)
        self.make_output_dirs()

        log('Loading selected gids')
        df_selected = self.input()['selected'].read() if not DEBUG else self.input()['selected'].read_cached()
        genera = sorted(set(df_selected.index))

        log(f'Collecting {len(genera):,} genera')
        if not DEBUG:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(worker, genera), total=len(genera)))
        else:
            results = [worker(x) for x in genera]
        df = pd.concat(results, ignore_index=True)

        if not DEBUG:
            self.save_hdf(df)


def worker(genus):
    path = os.path.join(DIR_OUT_FASTANI_INTER, 'for_random_non_reps', f'fastani_interspecies_to_{genus}.h5')
    df = pd.read_hdf(path)
    df['genus'] = genus
    return df
