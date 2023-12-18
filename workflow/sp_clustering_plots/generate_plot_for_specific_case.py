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

from workflow.external.gtdb_metadata import GtdbMetadataR207Full
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc
from workflow.util.log import log
from workflow.util.paths import get_gid_root

N_CPUS = mp.cpu_count()


def fastani_worker(gid, closest_rep_set, tmp_dir, gunc_source, max_css, domain):
    # Read the FASTA file
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')
    path_results_tmp = '/tmp/results2.h5'

    if os.path.isfile(path_results_tmp):
        return pd.read_hdf(path_results_tmp)

    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    if gunc_source == 'gtdb':
        path_contig_assign = os.path.join(gid_root, 'gunc_r95/gunc_output', f'{gid}.contig_assignments.tsv')
        gunc_domain = domain
    elif gunc_source == 'progenomes':
        path_contig_assign = os.path.join(gid_root, 'gunc_pro/gunc_output', f'{gid}.contig_assignments.tsv')
        if domain == 'd__Bacteria':
            gunc_domain = '2 Bacteria'
        elif domain == 'd__Archaea':
            gunc_domain = '2157 Archaea'
        else:
            raise ValueError(f'Unknown domain: {domain}')
    else:
        raise Exception(f'Unknown gunc source: {gunc_source}')
    df_contig_assign = read_contig_assignments_tsv(path_contig_assign)

    # Determine the percentage values at which this genome can have contigs removed
    taxon, tax_level = get_taxonomy_by_majority_vote_gunc(df_contig_assign, max_css, gunc_domain)
    d_pct_to_contigs_to_remove = contigs_to_remove_from_gunc(d_fna, df_contig_assign, taxon, tax_level)

    log(f'Running FastANI on {len(d_pct_to_contigs_to_remove):,} pct values')

    # Convert the rep set into their paths
    d_ref_paths = {x: os.path.join(get_gid_root(x), f'{x}.fna') for x in closest_rep_set}

    # Iterate over each percentage change
    cur_results = list()
    for pct, contigs_to_remove in tqdm(d_pct_to_contigs_to_remove.items()):

        # Write the fasta file without those contigs to disk
        qry_path = os.path.join(tmp_dir, f'{gid}_{pct}.fna')
        with open(qry_path, 'w') as f:
            SeqIO.write([x[1] for x in d_fna.items() if x[0] not in contigs_to_remove], f, 'fasta')

        # Prepare for FastANI
        ani = fastani(query=qry_path,
                      reference=list(d_ref_paths.values()),
                      cpus=N_CPUS,
                      single_execution=False,
                      bidirectional=True,
                      show_progress=False,
                      tmp_root=tempfile.gettempdir())

        # Prepare the ANI results for output
        d_ani = ani.as_dict()
        for ref_gid in d_ref_paths.keys():
            q_key = qry_path
            r_key = d_ref_paths[ref_gid]
            qvr = d_ani[q_key][r_key]
            rvq = d_ani[r_key][q_key]

            if qvr is not None and rvq is not None:
                ani = max(qvr.ani, rvq.ani)
                af = max(qvr.align_frac, rvq.align_frac)
            elif qvr is not None and rvq is None:
                ani = qvr.ani
                af = qvr.align_frac
            elif qvr is None and rvq is not None:
                ani = rvq.ani
                af = rvq.align_frac
            else:
                ani = 0
                af = 0
            cur_results.append((ref_gid, pct, ani, af))

    # Write the results to disk
    columns = ['reference', 'pct', 'ani', 'af']
    df = pd.DataFrame(cur_results, columns=columns)
    df.to_hdf(path_results_tmp, key='root', format='table', complevel=9, complib='blosc:lz4hc')
    return df


def main():

    # Klebsiella plot
    query_gid = 'GCA_900759445.1'
    ref_gids = ['GCF_000735435.1', 'GCF_001598295.1', 'GCF_002075345.1']

    # Blautia plot
    query_gid = 'GCA_900751995.1'
    ref_gids = ['GCF_001487165.1', 'GCF_015557635.1', 'GCF_000169235.1']
    # TOOD: Remember to remove temporary values

    log('Loading Max CSS')
    df_max_css = AggregateMaxCssLevelMerged().output().read_cached()

    log('Loading GTDB metadata')
    df_meta = GtdbMetadataR207Full().output().read_cached()

    log('Loading species clusters')
    df_sp_clusters = GtdbSpClustersR207().output().read_cached()

    css_row = df_max_css.loc[query_gid]

    source = css_row['source']

    # TODO: Automate this part
    domain = 'd__Bacteria'
    max_css = css_row['taxonomic_level']

    log('Calculating ANI values')
    with tempfile.TemporaryDirectory() as tmp_dir:
        pct_ani = fastani_worker(query_gid, set(ref_gids), tmp_dir, source, max_css, domain)
    base_ani = pd.read_hdf(os.path.join(get_gid_root(query_gid), 'fastani_gunc_failed_pct_0.h5'))

    df = pd.concat([base_ani, pct_ani], ignore_index=True)

    PCT_VALUES = [0, 1, 5, 10, 15, 20, 30, 40, 50]
    gid = query_gid

    meta_row = df_meta.loc[gid]
    cur_sp_radius = df_sp_clusters.loc[meta_row['species'], 'ani_radius']

    is_rep = meta_row['gtdb_genome_representative'][3:] == gid

    df = df[df['reference'] != gid]
    df = df[(df['pct'].isin(PCT_VALUES)) & (df['af'] >= 0.5) & (df['reference'].isin(ref_gids))]

    if len(df) == 0:
        print(f'no info left for: {gid}')
        return

    df = df.sort_values(by=['reference', 'pct'])

    # reps_that_exceed_radius = frozenset(df[df['ani'] >= cur_sp_radius]['reference'])
    d_rep_to_values = defaultdict(list)
    d_rep_to_pct = defaultdict(list)
    for _, row in df.iterrows():
        cur_ani = row['ani']
        cur_pct = row['pct']
        d_rep_to_values[row['reference']].append(cur_ani)
        d_rep_to_pct[row['reference']].append(cur_pct)

    cmap = sns.color_palette('bright', n_colors=len(PCT_VALUES), as_cmap=False)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.rcParams.update({'font.size': 12})

    plt.rcParams['svg.fonttype'] = 'none'
    next_colour_i = 0
    for i, (rep, ani) in enumerate(d_rep_to_values.items()):
        pct = d_rep_to_pct[rep]

        if df_meta.loc[rep, 'species'] == meta_row['species']:
            marker = 'o'
        else:
            marker = 'D'
        ax.plot(pct, ani, '-', marker=marker, label=df_meta.loc[rep, 'species'],
                color=cmap[next_colour_i % len(cmap)], zorder=20)
        next_colour_i += 1

        print(pct)
        print(ani)
        print(f'----------^^^  {df_meta.loc[rep, "species"]} ^^^-----------')

    ax.hlines(y=cur_sp_radius, xmin=-10, xmax=70, color='r',
              linestyle='--', alpha=0.4, label=f'ANI Radius: {cur_sp_radius:.2f}%')

    new_xticks, new_xlabels = list(), list()
    for i in PCT_VALUES:
        new_xticks.append(i)
        new_xlabels.append(str(i))

    new_yticks, new_ylabels = list(), list()
    for i in [-10] + PCT_VALUES:
        new_yticks.append(i)
        new_ylabels.append(str(i) if i not in {-10, 0} else '')

    ax.grid(True)
    ax.set_xlim(-1, 51)
    ax.set_ylim([94, 100])

    ax.set_xticklabels(new_xlabels)
    ax.set_xticks(new_xticks)
    ax.set_xlabel('% of genome removed')

    ax.set_ylabel('% ANI')

    if is_rep:
        title_text = f'{gid} (species representative)\n{meta_row["species"]}'
    else:
        title_text = f'{gid}\n{meta_row["species"]}'

    plt.title(title_text)

    plt.legend()

    ncbi_cat = df_meta.loc[gid, 'ncbi_genome_category']


    # os.makedirs('/tmp/guncplots', exist_ok=True)
    # plt.savefig(f"/tmp/guncplots/{gid}_rep_{is_rep}_{ncbi_cat}.svg")
    # plt.close()

    plt.savefig('/tmp/blautia.svg')
    # plt.show()

    print(gid)
    print(ncbi_cat)
    print('^^^^^^')
    return

    pass


if __name__ == '__main__':
    main()
