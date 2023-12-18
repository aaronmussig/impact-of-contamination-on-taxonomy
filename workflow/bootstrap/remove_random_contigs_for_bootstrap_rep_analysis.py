import multiprocessing as mp
import os
import tempfile
from collections import defaultdict, Counter

import pandas as pd
from magna.util.disk import move_file
from tqdm import tqdm

from workflow.bootstrap.remove_random_contigs_for_bootstrap_rep import RemoveRandomContigsForBootstrapRep
from workflow.bootstrap.remove_random_contigs_for_bootstrap_rep_sp_cluster import \
    RemoveRandomContigsForBootstrapRepSpCluster
from workflow.bootstrap.select_genomes_for_bootstrapping import SelectGenomesForBootstrapping
from workflow.config import DEBUG, DIR_OUT_BOOTSTRAP, PCT_VALUES
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
import matplotlib.pyplot as plt
import seaborn as sns



def worker(job):
    gid, repeat, gid_i, gid_rep, gid_species, ani_radius = job

    path_out = os.path.join(DIR_OUT_BOOTSTRAP, 'random_contigs_sp_clustering', str(repeat), f'{gid}_{gid_i}.h5')
    if os.path.isfile(path_out):
        os.remove(path_out)
        # return pd.read_hdf(path_out)

    path_fastani = os.path.join(DIR_OUT_BOOTSTRAP, 'random_contigs', str(repeat), f'{gid}_{gid_i}.h5')
    if not os.path.isfile(path_fastani):
        raise Exception(f'{path_fastani} does not exist ({gid} {gid_i} {repeat})')
    try:
        df_fastani = pd.read_hdf(path_fastani)
    except ValueError as e:
        if os.path.getsize(path_fastani) != 1024:
            raise Exception(f'Error parsing: {path_fastani}')
        # This is an empty file, likely no contigs were removed - skip it
        return None

    pct_values = frozenset(int(x) for x in df_fastani['pct'].unique())

    out = list()
    for pct_value in sorted(pct_values):
        df_fastani_subset = df_fastani[df_fastani['pct'] == pct_value]
        result = run_sp_clustering_on_fastani_values(gid, gid_rep, df_fastani_subset, ani_radius)
        result['pct'] = pct_value
        result['repeat'] = repeat
        result['gid_i'] = gid_i
        out.append(result)

    # Save intermediate file to disk
    os.makedirs(os.path.dirname(path_out), exist_ok=True)
    df = pd.DataFrame(out)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, 'out.h5')
        df.to_hdf(tmp_path, key='root', format='table', complevel=9, complib='blosc:lz4hc')
        move_file(tmp_path, path_out, checksum=True)
    return df


def calculate_bootstrap_statistic(df_meta, df_sp_clustering):

    # Some were originally run with more pct values, but just reduce to the analysis set
    df_changed = df_sp_clustering[df_sp_clustering['same'] == False]
    df_changed = df_changed[df_changed['pct'].isin(PCT_VALUES)]

    d_pct_to_changes = defaultdict(list)

    for cur_pct in PCT_VALUES:
        df_subset = df_changed[df_changed['pct'] == cur_pct]
        if len(df_subset) == 0:
            print(f'No values for pct {cur_pct}')
            continue

        df_subset = df_subset.sort_values(by=['gid', 'pct'], ascending=[True, False])

        for _, row in df_subset.iterrows():
            gid = row['gid']

            meta_row = df_meta.loc[gid]
            is_rep = meta_row['gtdb_representative'] == 't'

            if is_rep:
                change_type = 'representative changed species'
            else:
                if row['type'] == 'no_af' or row['type'] == 'no_ani':
                    change_type = 'novel_sp_cluster'
                else:
                    change_type = 'species changed'

            d_pct_to_changes[cur_pct].append(change_type)

    return d_pct_to_changes



class RemoveRandomContigsForBootstrapRepAnalysis(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'sp_clustering': RemoveRandomContigsForBootstrapRepSpCluster(),
            'sp_clusters': GtdbSpClustersR207(),
            'meta': GtdbMetadataR207(),
            'bootstrap_gids': SelectGenomesForBootstrapping(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_BOOTSTRAP, 'remove_random_contigs_for_bootstrap_plot_data.h5'))

    def run(self):
        log('Running species clustering on gunc failed contig by bootstrap random analysis', title=True)
        # self.make_output_dirs()

        log('Loading bootstrap genomes')
        df_bootstrap = self.input()['bootstrap_gids'].read() if not DEBUG else self.input()['bootstrap_gids'].read_cached()

        log('Loading GTDB metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        # log('Loading GTDB species clusters')
        # df_sp_clusters = self.input()['sp_clusters'].read() if not DEBUG else self.input()['sp_clusters'].read_cached()

        log('Loading species clustering results')
        df_sp_clustering = self.input()['sp_clustering'].read() if not DEBUG else self.input()['sp_clustering'].read_cached()

        # Calculate the statistic for each bootstrap repetition
        d_repeat_to_values = dict()
        for bootstrap_repeat, genomes in df_bootstrap.iterrows():
            df_subset = df_sp_clustering[df_sp_clustering['repeat'] == bootstrap_repeat]
            result = calculate_bootstrap_statistic(df_meta, df_subset)
            d_repeat_to_values[bootstrap_repeat] = result

        # Average the results out over each bootstrap repetition
        plot_data = plot_results(d_repeat_to_values)

        self.save_hdf(plot_data)

        return



def plot_results(d_repeat_to_values, n_repeats=10):

    # Convert them to counts
    d_repeat_to_pct_changes = defaultdict(dict)
    for repeat, d_pct_to_changes in d_repeat_to_values.items():
        for pct, changes in d_pct_to_changes.items():
            d_repeat_to_pct_changes[repeat][pct] = Counter(changes)

    # Get the average value to plot as a bar graph
    d_pct_to_counts = defaultdict(lambda: defaultdict(lambda: 0))
    for repeat, d_pct_to_changes in d_repeat_to_pct_changes.items():
        for pct, d_change_to_count in d_pct_to_changes.items():
            for change, count in d_change_to_count.items():
                d_pct_to_counts[pct][change] += count

    # Average the results over each repetition
    d_pct_to_counts_avg = defaultdict(dict)
    for repeat, d_pct_to_changes in d_repeat_to_pct_changes.items():
        for pct, d_change_to_count in d_pct_to_changes.items():
            for change, count in d_change_to_count.items():
                d_pct_to_counts_avg[pct][change] = count / n_repeats

    # Generate the data
    rows = list()
    for pct, d_change_to_count in d_pct_to_counts_avg.items():
        for change, count in d_change_to_count.items():
            rows.append({'pct': pct, 'change': change, 'count': count})
    df_avg = pd.DataFrame(rows)

    d_change_to_text = {
        'species changed': 'Genome changed species',
        'representative changed species': 'Species clusters merged',
        'novel_sp_cluster': 'New species cluster formed',
    }

    df_avg['change'] = df_avg['change'].apply(lambda x: d_change_to_text[x])

    df_avg.sort_values(by='change', ascending=False, inplace=True)

    fig, ax = plt.subplots(figsize=(15, 10))
    plt.rcParams['svg.fonttype'] = 'none'
    ax.grid(True)
    plt.rcParams.update({'font.size': 18})
    sns.barplot(data=df_avg, x='pct', y='count', hue='change', ax=ax)
    ax.set_ylabel('Number of putatively contaminated genomes affected')
    ax.set_xlabel('% of genome removed')
    ax.set_ylim(0, 12)
    ax.set_yticks(list(range(0, 12)))

    plt.title('Number of putatively contaminated genomes affected by percentage of genome removed')

    plt.savefig('/tmp/bootstrap.svg')
    plt.show()

    return df_avg


