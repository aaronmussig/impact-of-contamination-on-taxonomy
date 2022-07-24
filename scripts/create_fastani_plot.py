from collections import defaultdict

from workflow.external.gtdb_metadata import GtdbMetadataR207

from chiton.fastani import fastani
import numpy as np

from workflow.util.paths import get_gid_root
import seaborn as sns
import matplotlib.pyplot as plt

import os
import pandas as pd

PATH_TMP = '/tmp/fastani_results.tsv'
# sns.set_theme(style="viridis", palette=None)

def gen_data():
    rank = 'g__Helicobacter'

    df_meta = GtdbMetadataR207().output().read_cached()
    df_subset = df_meta[(df_meta['genus'] == rank) & (df_meta['gtdb_representative'] == 't')]

    gids = {'GCF_002026585.1'}

    gid_to_species = dict()
    for gid, row in df_subset.iterrows():
        gid_to_species[gid] = row['species']

    gid_to_species['GCF_002026585.1'] = 's__Helicobacter pylori_CI'

    gid_to_path = dict()
    path_to_gid = dict()
    for gid, _ in gid_to_species.items():
        gid_to_path[gid] = os.path.join(get_gid_root(gid), f'{gid}.fna')
        path_to_gid[gid_to_path[gid]] = gid

    # Calculate ANI
    ani = fastani(query=list(gid_to_path.values()),
                  reference=list(gid_to_path.values()),
                  cpus=10,
                  single_execution=False,
                  bidirectional=True)

    # Prepare the ANI results for output
    d_ani = ani.as_dict()
    rows = list()
    for qry_path, d_vals in d_ani.items():
        for ref_path, result in d_vals.items():
            rows.append({
                'query': path_to_gid[qry_path],
                'reference': path_to_gid[ref_path],
                'ani': result.ani,
                'af': result.align_frac
            })

    df = pd.DataFrame(rows)

    df.to_csv(PATH_TMP, sep='\t', index=False)
    return



def main():

    if not os.path.isfile(PATH_TMP):
        gen_data()

    df_meta = GtdbMetadataR207().output().read_cached()

    df = pd.read_csv(PATH_TMP, sep='\t')

    d_results = defaultdict(dict)
    for _, row in df.iterrows():
        d_results[row['query']][row['reference']] = row['ani']

    # Create a matrix
    keys = sorted(d_results.keys())
    mat = np.zeros((len(keys), len(keys)))
    for i, key1 in enumerate(keys):
        for j, key2 in enumerate(keys):
            mat[i, j] = d_results[key1][key2]
    for i, key1 in enumerate(keys):
        for j, key2 in enumerate(keys):
            mat[i, j] = max(mat[i, j], mat[j, i])

    species = [df_meta.loc[x, 'species'] if x != 'GCF_002026585.1' else 's__Helicobacter pylori_CI *' for x in keys]

    df = pd.DataFrame(mat, columns=species)

    df['idx'] = species
    df = df.set_index('idx')


    plt.figure(figsize=(10, 6))
    plt.rcParams['svg.fonttype'] = 'none'
    g = sns.clustermap(df, annot=True, cmap='viridis',  fmt='.1f')

    # mask = np.tril(np.ones_like(df))
    # values = g.ax_heatmap.collections[0].get_array().reshape(df.shape)
    # new_values = np.ma.array(values, mask=mask)
    # g.ax_heatmap.collections[0].set_array(new_values)

    plt.savefig('/tmp/heatmap.svg')
    # plt.show()

    return




if __name__ == '__main__':
    main()