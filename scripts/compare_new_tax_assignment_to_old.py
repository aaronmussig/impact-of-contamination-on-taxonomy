from magna.gunc import read_contig_assignments_tsv
from tqdm import tqdm

from workflow.config import DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.get_taxonomy_from_gunc_contig_assignment import GetTaxonomyFromGuncContigAssignment
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc
from workflow.util.paths import get_gid_root
import multiprocessing as mp

import os

def worker(job):
    gid, max_css, expected_domain = job

    gid_root = get_gid_root(gid)

    gunc_dir_r95 = os.path.join(gid_root, 'gunc_r95')
    gunc_dir_pro = os.path.join(gid_root, 'gunc_pro')

    out = list()
    for gunc_source, gunc_root in [('progenomes', gunc_dir_pro), ('gtdb', gunc_dir_r95)]:
        df_contig_assign = read_contig_assignments_tsv(os.path.join(gunc_root, 'gunc_output', f'{gid}.contig_assignments.tsv'))

        if gunc_source == 'progenomes':
            if expected_domain == 'd__Archaea':
                domain = '2157 Archaea'
            else:
                domain = '2 Bacteria'
        else:
            domain = expected_domain

        taxon, tax_level = get_taxonomy_by_majority_vote_gunc(df_contig_assign, max_css, domain)
        out.append({
            'gid': gid,
            'db': gunc_source,
            'taxon': taxon,
            'tax_level': tax_level,
        })

    return out


def main():
    df = GetTaxonomyFromGuncContigAssignment().output().read_cached()

    df_meta = GtdbMetadataR207().output().read_cached()
    df_max_css = AggregateMaxCssLevelMerged().output().read_cached()

    queue = list()
    for gid, row in df_max_css.iterrows():
        max_css = row['taxonomic_level']
        expected_domain = df_meta.loc[gid]['domain']
        queue.append((gid, max_css, expected_domain))

    if DEBUG:
        queue = queue[0:10]

    rows = list()
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))
    for result in results:
        rows.extend(result)

    ranks = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    for row in rows:
        df_tax = df[(df['gid'] == row['gid']) & (df['db'] == row['db'])]['taxonomy'].values[0]

        highest_level_new = ranks[len(df_tax.split(';'))]
        if not df_tax.endswith(row['taxon']):
            print(row['gid'], row['db'], row['taxon'], df_tax)

    return



if __name__ == '__main__':
    main()