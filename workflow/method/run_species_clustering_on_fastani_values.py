import os
from typing import Optional

import pandas as pd

from workflow.config import R207_AF_THRESHOLD
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def run_sp_clustering_on_fastani_values(
        gid: str,
        rep_gid: str,
        df_fastani: pd.DataFrame,
        ani_radius: float,
        row_id: Optional[str] = None,
):
    # Alternative row id for returning the genome id
    if row_id is None:
        row_id = gid

    # If this is a species representative, then remove hits to self
    is_rep = gid == rep_gid
    if is_rep:
        # log(f'Removing self hits for {gid} as it is a species representative')
        df_fastani = df_fastani[df_fastani['reference'] != gid]

    # Filter to only those within the AF threshold
    df_fastani = df_fastani[df_fastani['af'] >= R207_AF_THRESHOLD]

    # No genomes are within AF threshold = novel species cluster
    if len(df_fastani) == 0:
        # log(f'No genomes within AF radius {gid}')
        return {
            'gid': row_id,
            'new_sp_rep': None,
            'ani': None,
            'af': None,
            'type': 'no_af',
            'same': is_rep
        }

    # Filter by ANI radius
    df_fastani = df_fastani[df_fastani['ani'] >= ani_radius]

    # No genomes within ANI radius = novel species cluster
    if len(df_fastani) == 0:
        # log(f'No genomes within ANI radius {gid}')
        return {
            'gid': row_id,
            'new_sp_rep': None,
            'ani': None,
            'af': None,
            'type': 'no_ani',
            'same': is_rep
        }

    # Successfully found a species representative
    # Take the most significant hit
    df_fastani = df_fastani.sort_values(by=['ani', 'af'], ascending=[False, False])
    new_sp_rep_series = df_fastani.iloc[0]
    new_sp_rep_gid = new_sp_rep_series['reference']

    # if new_sp_rep['reference'] in failed_gids:
    #     raise NotImplemented(f'New representative is chimeric: {gid} {cutoff}')

    if new_sp_rep_gid != rep_gid:
        # log(f'new sp rep is not the same as old sp rep: {gid}')
        is_same = False
    else:
        is_same = True

    return {
        'gid': row_id,
        'new_sp_rep': new_sp_rep_gid,
        'ani': float(new_sp_rep_series['ani']),
        'af': float(new_sp_rep_series['af']),
        'type': 'sp_rep',
        'same': is_same
    }


def _main():
    gid = 'GCA_014894375.1'

    gid_root = get_gid_root(gid)
    df_fastani = pd.read_hdf(os.path.join(gid_root, 'fastani_gunc_failed_by_contamination.h5'))

    df_sp_clusters = GtdbSpClustersR207().output().read_cached()

    df_meta = GtdbMetadataR207().output().read_cached()

    gid_rep = df_meta.loc[gid, 'gtdb_genome_representative'][3:]
    gid_species = df_meta.loc[gid, 'species']

    pct_values = frozenset(int(x) for x in df_fastani['pct'].unique())

    ani_radius = float(df_sp_clusters.loc[gid_species, 'ani_radius'])

    for pct_value in pct_values:
        df_fastani_subset = df_fastani[df_fastani['pct'] == pct_value]
        result = run_sp_clustering_on_fastani_values(gid, gid_rep, df_fastani_subset, ani_radius)
        print(result)

    return


if __name__ == '__main__':
    _main()
