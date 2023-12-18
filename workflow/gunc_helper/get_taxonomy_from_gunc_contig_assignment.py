import multiprocessing as mp
import os
from collections import Counter

import pandas as pd
from gunc.get_scores import create_base_data as gunc_create_base_data
from gunc.get_scores import get_abundant_lineages_cutoff as gunc_get_abundant_lineages_cutoff
from gunc.get_scores import get_stats as gunc_get_stats
from gunc.get_scores import read_diamond_output as gunc_read_diamond_output
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_GUNC
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root

RANKS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


def worker(job):
    gid, max_css, source_db = job
    gid_root = get_gid_root(gid)

    if source_db == 'gtdb':
        gunc_dir = os.path.join(gid_root, 'gunc_r95')
        gunc_db_str = 'gtdb_95'
    elif source_db == 'progenomes':
        gunc_dir = os.path.join(gid_root, 'gunc_pro')
        gunc_db_str = 'progenomes_2.1'
    else:
        raise ValueError(f'Unknown source_db: {source_db}')

    # This repeats the method used in the GUNC program
    diamond_file_path = os.path.join(gunc_dir, 'diamond_output', f'{gid}.diamond.{gunc_db_str}.out')
    diamond_df = gunc_read_diamond_output(diamond_file_path)
    base_data = gunc_create_base_data(diamond_df, gunc_db_str)

    genes_mapped, contig_count = gunc_get_stats(diamond_df)
    abundant_lineages_cutoff = gunc_get_abundant_lineages_cutoff(False, genes_mapped)

    tax_data = base_data.groupby(max_css).filter(
        lambda x: len(x) > abundant_lineages_cutoff
    )

    df_taxstrings = tax_data[list(RANKS[:RANKS.index(max_css) + 1])]
    taxstrings = [';'.join(x) for x in df_taxstrings.values]
    taxstrings_count = Counter(taxstrings)

    # Can't reach a consensus if there are multiple most common
    if len(set(x[1] for x in taxstrings_count.most_common(2))) == 1:
        raise Exception(f'{gid} has multiple most common taxstrings')

    # Get the most common taxstring
    taxonomy = taxstrings_count.most_common(1)[0][0]

    return {
        'gid': gid,
        'db': source_db,
        'taxonomy': taxonomy,
    }


class GetTaxonomyFromGuncContigAssignment(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_GUNC, 'taxonomy_from_contig_assignment.h5'))

    def run(self):
        log('Getting taxonomy from GUNC contig assignment', title=True)
        self.make_output_dirs()

        log('Loading MaxCSS file')
        df_max_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Creating queue')
        queue = list()
        for gid, row in tqdm(df_max_css.iterrows(), total=len(df_max_css)):
            # if gid != 'GCF_001037485.2':
            #     continue

            queue.append((gid, row['taxonomic_level'], row['source']))

        log('Processing queue')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        df = pd.DataFrame(results)
        if not DEBUG:
            self.save_hdf(df, index='gid')
        return
