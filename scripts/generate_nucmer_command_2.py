from collections import defaultdict
from typing import List

from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.util.paths import get_gid_root, get_gid_r207_root, get_gid_r207_fna
import os
import pandas as pd


GIDS = ['GCA_900751995.1', 'GCA_900759445.1', 'GCF_002026585.1', 'GCA_003086435.1', 'GCA_900555355.1']

def run_nucmer(ref_gid: str, qry_gids: List[str]):


    gid_fna =get_gid_r207_fna(ref_gid)

    cmd = [
        'mummer2circos',
        '-l',
        '-r',
        gid_fna,
        '-q',
    ]

    for qry_gid in qry_gids:
        qry_gid_fna =get_gid_r207_fna(qry_gid)
        cmd.append(qry_gid_fna)

    print(' '.join(cmd))
    return


def main():

    run_nucmer('GCA_016712675.1', qry_gids=[
        'GCA_016703265.1', # from genus LZORAL
        'GCF_012939995.1' # from competibacter
    ])

    print('')

    # for contaminatn

    # run_nucmer('GCF_900291955.1', qry_gids=['GCF_001487165.1', 'GCA_014385135.2'])
    pass





if __name__ == '__main__':
    main()