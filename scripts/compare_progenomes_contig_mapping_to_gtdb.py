import multiprocessing as mp
import os

import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, PCT_VALUES
from workflow.cutoff.create_cutoff_from_max_css import CreateCutoffPctFromMaxCss
from workflow.external.gtdb_metadata import GtdbMetadataR95
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR95
from workflow.fastani.run_species_clustering_on_failed_aggregate import RunSpeciesClusteringOnFailedAggregate
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.util.paths import get_gid_root_pct


if DEBUG:
    DF = RunSpeciesClusteringOnFailedAggregate().output().read_cached()
    DF_META = GtdbMetadataR95().output().read_cached()
    DF_MAX_CSS = AggregateMaxCssLevelMerged().output().read_cached()
    DF_CUTOFF = CreateCutoffPctFromMaxCss().output().read_cached()
    DF_SP_CLUSTERS = GtdbSpClustersR95().output().read_cached()
else:
    DF = RunSpeciesClusteringOnFailedAggregate().output().read()
    DF_META = GtdbMetadataR95().output().read()
    DF_MAX_CSS = AggregateMaxCssLevelMerged().output().read()
    DF_CUTOFF = CreateCutoffPctFromMaxCss().output().read()
    DF_SP_CLUSTERS = GtdbSpClustersR95().output().read()

def worker(gid):
    gid_sp_clustering_df = DF[DF['gid'] == gid]
    meta_row = DF_META.loc[gid]
    is_rep = meta_row['gtdb_representative'] == 't'
    is_type = meta_row['gtdb_type_designation'] != 'not type material'

    d_pct_to_changes = dict()
    df_changed = gid_sp_clustering_df[gid_sp_clustering_df['same'] == False]

    if len(df_changed) == 0:
        return

    for pct in PCT_VALUES:
        df_subset = df_changed[df_changed['cutoff'] <= pct]
        if len(df_subset) == 0:
            continue

        for _, row in df_subset.iterrows():
            if is_rep:
                if not is_type:
                    d_pct_to_changes[pct] = ('representative changed species', row)
                else:
                    if row['type'] == 'no_af' or row['type'] == 'no_ani':
                        d_pct_to_changes[pct] = ('novel_sp_cluster', row)
                    else:
                        d_pct_to_changes[pct] = ('species changed', row)

    if len(d_pct_to_changes) == 0:
        return

    css_row = DF_MAX_CSS.loc[gid]
    gunc_source = css_row['source']

    cutoff_row = DF_CUTOFF.loc[gid]

    out = list()
    for pct, (change, row) in sorted(d_pct_to_changes.items(), key=lambda x: x[0]):
        out.append({
            'gid': gid,
            'gunc_db': gunc_source,
            'pct_genome_removed': pct,
            'type_of_change': change,
            'expected_species': meta_row['species'],
            'new_species': row['new_sp_rep'],
            'ani': row['ani'],
            'af': row['af'],
            'is_identical': row['same']
        })
    return out


def worker2(gid):
    # gid_sp_clustering_df = DF[DF['gid'] == gid]
    meta_row = DF_META.loc[gid]
    is_rep = meta_row['gtdb_representative'] == 't'
    is_type = meta_row['gtdb_type_designation'] != 'not type material'
    # cutoff_row = DF_CUTOFF.loc[gid]
    cur_rep_gid = meta_row['gtdb_genome_representative'][3:]
    css_row = DF_MAX_CSS.loc[gid]
    gunc_source = css_row['source']
    ani_df = DF[(DF['gid'] == gid) & (DF['cutoff'] > 0)]

    changed_at_any_pct = len(ani_df[ani_df['same'] == False]) > 0

    cur_sp_cluster = DF_SP_CLUSTERS.loc[meta_row['species']]


    out = list()
    for _, row in ani_df.iterrows():
        cur_pct = int(row['cutoff'])

        cur_ani_df_row = ani_df[ani_df['cutoff'] == cur_pct].iloc[0]

        gid_root_pct = get_gid_root_pct(gid, cur_pct)
        cur_fastani_path = os.path.join(gid_root_pct, 'fastani', 'closest_rep_ani_af.tsv')
        cur_fastani_raw_df = pd.read_csv(cur_fastani_path, sep='\t')

        cur_rep_hit = cur_fastani_raw_df[cur_fastani_raw_df['reference'] == cur_rep_gid].iloc[0]

        if cur_ani_df_row['same'] == True:
            type_of_change = 'none'
        else:
            if cur_ani_df_row['type'] == 'sp_rep':
                type_of_change = 'changed species cluster'
            elif cur_ani_df_row['type'] in {'no_af', 'no_ani'}:
                type_of_change = 'new species cluster'
            else:
                raise ValueError('Unknown type: {}'.format(cur_ani_df_row['type']))

        out.append({
            'gid': gid,
            'changed_at_any_pct': changed_at_any_pct,
            'is_gtdb_rep': is_rep,
            'is_type_material': is_type,
            'pct_genome_removed': cur_pct,
            'gunc_db': gunc_source,

            'type_of_change': type_of_change,

            'expected_species': meta_row['species'],
            'gtdb_r95_ani_radius': cur_sp_cluster['ani_radius'],

            'new_sp_rep': cur_ani_df_row['new_sp_rep'],

            'ani_to_new_sp_rep': cur_ani_df_row['ani'],
            'af_to_new_sp_rep': cur_ani_df_row['af'],

            'ani_to_gtdb_sp_rep': cur_rep_hit['ani'],
            'af_to_gtdb_sp_rep': cur_rep_hit['af'],

        })

    return out


def main():
    gids = sorted(set(DF_MAX_CSS.index))

    if DEBUG:
        results = [worker2(x) for x in gids if x == 'GCF_000421645.1']
    else:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(worker2, gids), total=len(gids)))

    rows = list()
    [rows.extend(r) for r in results if r is not None]

    df = pd.DataFrame(rows)

    df.sort_values(by=['gid', 'pct_genome_removed'], inplace=True)
    df.to_csv('/srv/home/uqamussi/sp_clustering.tsv', sep='\t', index=False)

    pass


if __name__ == '__main__':
    main()
