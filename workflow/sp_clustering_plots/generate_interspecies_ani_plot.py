import multiprocessing as mp
import os
import tempfile
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from chiton.fastani import fastani
from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.external.gtdb_metadata import GtdbMetadataR207Full, GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc
from workflow.util.log import log
from workflow.util.paths import get_gid_root



def main():
    log('Loading metadata')
    df_meta = GtdbMetadataR207().output().read_cached()
    df_meta = df_meta[df_meta['gtdb_representative'] == 't']

    log('Mapping genus to species representatives')
    d_genus_to_species_reps = defaultdict(set)
    for gid, row in tqdm(df_meta.iterrows(), total=len(df_meta)):
        d_genus_to_species_reps[row['genus']].add(gid)

if __name__ == '__main__':
    main()
