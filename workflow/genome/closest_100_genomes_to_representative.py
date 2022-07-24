import os
from typing import Dict, Tuple

import dendropy
import numpy as np
import pandas as pd
from phylodm import PhyloDM
from tqdm import tqdm

from workflow.config import DIR_OUT_GENOMES
from workflow.external.gtdb_tree import GtdbTreeArcR207, GtdbTreeBacR207
from workflow.genome.symlink import GenomeSymlink
from workflow.model.luigi import LocalTargetHdf5
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


def get_closest_n_reps_for_gids(n: int, path_tree) -> Dict[str, Tuple[str]]:
    tree = dendropy.Tree.get(path=path_tree, schema='newick', preserve_underscores=True)

    print('Loading tree with PhyloDM')
    phylodm = PhyloDM.load_from_dendropy(tree)

    print('Computing distance matrix')
    dm = phylodm.dm()
    pdm_taxa = tuple(x[3:] for x in phylodm.taxa())

    print('Formatting representatives')
    d_sp_rep_to_closest = dict()
    for idx, sp_rep in tqdm(enumerate(pdm_taxa), total=len(pdm_taxa)):
        # Get the closest 100
        closest_n_idx = np.argsort(dm[idx])[:n]

        d_sp_rep_to_closest[sp_rep] = list()
        for closest_idx in closest_n_idx:
            d_sp_rep_to_closest[sp_rep].append(pdm_taxa[closest_idx])

        if len(d_sp_rep_to_closest[sp_rep]) != 100:
            raise ValueError('Expected 100 closest reps for {}'.format(sp_rep))

    return {k: tuple(v) for k, v in d_sp_rep_to_closest.items()}


class Closest100GenomesToRepresentative(LuigiTask):

    def requires(self):
        return {
            '_genomes': GenomeSymlink(),
            'ar122_tree': GtdbTreeArcR207(),
            'bac120_tree': GtdbTreeBacR207(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_GENOMES, 'closest_100_representatives.h5'))

    def run(self):
        log('Determining closest representative genomes', title=True)
        self.make_output_dirs()

        path_arc_tree = self.input()['ar122_tree'].path
        path_bac_tree = self.input()['bac120_tree'].path

        arc_sp_rep_closest = get_closest_n_reps_for_gids(100, path_arc_tree)
        bac_sp_rep_closest = get_closest_n_reps_for_gids(100, path_bac_tree)
        all_sp_rep_closest = {**arc_sp_rep_closest, **bac_sp_rep_closest}

        print(f'Creating dataframe rows for {len(all_sp_rep_closest):,} representatives')
        rows = list()
        for gid, lst_closest in all_sp_rep_closest.items():
            rows.append({'gid': gid, 'closest_representatives': '|'.join(lst_closest)})

        print('Creating dataframe')
        df = pd.DataFrame(rows)
        self.save_hdf(df, index='gid')
