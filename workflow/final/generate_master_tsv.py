import os
from collections import defaultdict

import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
from workflow.config import DIR_OUT_FINAL
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani_contig_split.b_run_fastani import FastAniContigSplitRunFastAni
from workflow.fastani_contig_split.d_report_results import FastAniContigSplitReportResultsFastAni
from workflow.fastani_contig_split.d_report_results_relaxed import FastAniContigSplitReportResultsFastAniRelaxed
from workflow.fasttree_marker_split.e_decorate_fasttree import FastTreeMarkerSplitDecorateFastTree
from workflow.fasttree_marker_split.f_a_compare_patristic_distance import FastTreeMarkerSplitComparePatristicDistance
from workflow.fasttree_marker_split.f_b_calculate_num_jumps_between_halves import FastTreeMarkerSplitJumpsBetweenHalves
from workflow.fasttree_marker_split_true_case.c_calculate_patristic_distance import \
    FastTreeMarkerSplitTrueCaseCalculatePatristicDistance
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetTsv
from workflow.util.log import log
from workflow.util.taxonomy import calculate_taxonomic_novelty, D_TAX_NOVELTY_TO_INDEX
from workflow.util.tree import parse_label
from scipy.stats import lognorm

RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')

RANKS_STRAIN = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')

N_MARKER_SPLIT_TRUE = 10

class FinalGenerateMasterTsv(LuigiTask):

    def requires(self):
        req =  {
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
            'marker_split': FastTreeMarkerSplitDecorateFastTree(),
            'marker_split_pd': FastTreeMarkerSplitComparePatristicDistance(),
            # 'fastani_split': FastAniContigSplitReportResultsFastAni(),
            'fastani_split': FastAniContigSplitReportResultsFastAniRelaxed(),
            'marker_split_jumps': FastTreeMarkerSplitJumpsBetweenHalves(),
            'ani_raw': FastAniContigSplitRunFastAni(),
        }
        for i in range(N_MARKER_SPLIT_TRUE):
            req[f'marker_split_true_pd_{i}'] = FastTreeMarkerSplitTrueCaseCalculatePatristicDistance(replicate_id=i)
        return req

    def output(self):
        return LocalTargetTsv(os.path.join(DIR_OUT_FINAL, 'master.tsv'))

    def load_marker_split_taxonomy(self):
        out = dict()
        path = self.input()['marker_split'].path
        path = os.path.join(os.path.dirname(path), 'decorated.tree-taxonomy')
        with open(path) as f:
            for line in f.readlines():
                gid, tax = line.strip().split('\t')
                out[gid] = tax
        return out

    def load_marker_split_true_pd(self):

        rows = list()

        for i in range(N_MARKER_SPLIT_TRUE):
            df = self.input()[f'marker_split_true_pd_{i}'].maybe_read_cached()

            for row in df.itertuples():
                rows.append({
                    'replicate_id': i,
                    'novelty': row.tax_novelty,
                    'pd': row.pd,
                })

        df_out = pd.DataFrame(rows)
        df_out['pd'] = df_out['pd'].apply(lambda x: x * 100000000)






        plot = False

        if plot:
            import seaborn as sns
            import matplotlib.pyplot as plt
            sns.boxplot(df_out, x='novelty', y='pd', hue='replicate_id')
            plt.show()

        strain_loc, strain_scale, strain_mean, strain_std = calculate_for_strain(df_out, 'strain', 1, 10000, 0, 10, plot=plot)
        species_loc, species_scale, species_mean, species_std = calculate_for_strain(df_out, 'species', 1, plot_min=0, plot_max=15, plot=plot)
        genus_loc, genus_scale, genus_mean, genus_std = calculate_for_strain(df_out, 'genus', plot_min=0, plot_max=15, plot=plot)


        return {
            'strain': {
                'loc': strain_loc,
                'scale': strain_scale,
                'mean': strain_mean,
                'std': strain_std,
            },
            'species': {
                'loc': species_loc,
                'scale': species_scale,
                'mean': species_mean,
                'std': species_std,
            },
            'genus': {
                'loc': genus_loc,
                'scale': genus_scale,
                'mean': genus_mean,
                'std': genus_std,
            },
        }

    def load_top_hit_from_ani_raw_data(self):
        df_ani_raw = self.input()['ani_raw'].maybe_read_cached()

        out = dict()
        for row in tqdm(df_ani_raw.itertuples(), total=len(df_ani_raw)):
            gid = row.query
            if gid.replace('_C', '') == row.ref:
                continue
            if gid not in out:
                out[gid] = {
                    'ref': row.ref,
                    'ani': row.ani,
                    'af': row.af
                }
            else:
                if row.ani > out[gid]['ani']:
                    out[gid] = {
                        'ref': row.ref,
                        'ani': row.ani,
                        'af': row.af
                    }
        return out


    def run(self):
        log(f'FinalGenerateMasterTsv', title=True)
        self.make_output_dirs()

        log('Collecting raw FastANI data')
        d_gid_to_ani_raw_top_hit = self.load_top_hit_from_ani_raw_data()

        log('Loading marker split true pd')
        d_distributions = self.load_marker_split_true_pd()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()
        d_gid_to_contig_count = df_meta['contig_count'].to_dict()
        d_tax_to_gids = defaultdict(set)
        for gid, tax in tqdm(d_gid_to_tax.items()):
            d_tax_to_gids[tax].add(gid)

        log('Loading MaxCSS')
        df_css = self.input()['max_css'].maybe_read_cached()
        fail_gids = {x for x in df_css.index if d_gid_to_tax[x].startswith('d__Bacteria')}

        d_gid_to_sp_where_all_fail = get_species_where_all_fail(df_css, df_meta)

        log('Loading marker split taxonomy')
        d_gid_to_marker_split_tax = self.load_marker_split_taxonomy()

        log('Loading fastani split results')
        df_fastani = self.input()['fastani_split'].maybe_read_cached()
        df_fastani.set_index('gid', inplace=True)

        log('Getting taxonomic novelty (subset to the genomes in the tree)')
        gids_in_tree = {x.replace('TEST_', '').replace('_C', '') for x in d_gid_to_marker_split_tax}
        d_gid_to_tax_novelty, d_tax_novelty_to_gid = calculate_taxonomic_novelty(
            {x: d_gid_to_tax[x] for x in gids_in_tree}
        )

        log('Loading patristic distances (for marker split tree)')
        df_pd = self.input()['marker_split_pd'].maybe_read_cached()
        df_pd['pd_between_halves'] = df_pd['pd_between_halves'].apply(lambda x: x * 100000000)
        df_pd.set_index('gid', inplace=True)
        d_gid_to_pd = df_pd['pd_between_halves'].to_dict()

        log('Loading marker split tree')
        tree = self.input()['marker_split'].read()

        # log('Encoding tree bipartitions')
        # tree.encode_bipartitions()

        log('Calculating number of nodes between chimeric halves')
        df_marker_split_jumps = self.input()['marker_split_jumps'].maybe_read_cached()

        log('Reading taxonomy from tree')
        d_gid_to_tree_tax = get_taxonomy_from_tree(tree)

        log('Calculating the highest rank that agrees for each failed genome')
        d_gid_to_highest_agree = calculate_highest_agree_rank(
            fail_gids,
            d_gid_to_tax_novelty,
            d_gid_to_tax,
            d_gid_to_tree_tax
        )

        log('Creating output dataframe')
        rows = list()

        for gid in tqdm(sorted(fail_gids)):

            keep_highest_agree = d_gid_to_highest_agree[gid]
            chim_highest_agree = d_gid_to_highest_agree[f'{gid}_C']

            gid_keep_split_tax_parsed = d_gid_to_tree_tax[gid]
            gid_chim_split_tax_parsed = d_gid_to_tree_tax[f'{gid}_C']

            lst_tax_keep_half = list()
            for rank in RANKS[RANKS.index(keep_highest_agree):]:
                lst_tax_keep_half.extend(gid_keep_split_tax_parsed[rank[0]])

            lst_tax_chim_half = list()
            for rank in RANKS[RANKS.index(chim_highest_agree):]:
                lst_tax_chim_half.extend(gid_chim_split_tax_parsed[rank[0]])

            novelty = d_gid_to_tax_novelty[gid]

            if novelty == 'strain':
                keep_is_chim = keep_highest_agree != 'species'
                chim_is_chim = chim_highest_agree != 'species'
            else:
                novelty_idx = RANKS.index(novelty)
                keep_is_chim = RANKS.index(keep_highest_agree) < novelty_idx
                chim_is_chim = RANKS.index(chim_highest_agree) < novelty_idx

            fastani_row = df_fastani.loc[gid]

            ani_keep_tax_lst = list()
            if str(fastani_row['keep_tax']) != 'nan':
                for true_rank, test_rank in zip(d_gid_to_tax[gid].split(';'), fastani_row['keep_tax'].split(';')):
                    if true_rank != test_rank:
                        ani_keep_tax_lst.append(test_rank)

            ani_chim_tax_lst = list()
            if str(fastani_row['disc_tax']) != 'nan':
                for true_rank, test_rank in zip(d_gid_to_tax[gid].split(';'), fastani_row['disc_tax'].split(';')):
                    if true_rank != test_rank:
                        ani_chim_tax_lst.append(test_rank)

            if fastani_row['keep_type'] == 'sp_rep':
                ani_keep_half_is_wrong = 'no' if fastani_row['keep_same_as_207'] else 'yes'
            else:
                ani_keep_half_is_wrong = 'inconclusive'

            if fastani_row['disc_type'] == 'sp_rep':
                ani_chim_half_is_wrong = 'no' if fastani_row['disc_same_as_207'] else 'yes'
            else:
                ani_chim_half_is_wrong = 'inconclusive'

            n_jumps = df_marker_split_jumps.loc[gid, 'n_jumps']

            gunc_row = df_css.loc[gid]

            # Calculate how many standard deviations this point is away from the mean
            if novelty in {'strain', 'species', 'genus'}:
                n_std_away = str(round((np.log(d_gid_to_pd[gid]) - d_distributions[novelty]['mean']) / d_distributions[novelty]['std'], 4))
            else:
                n_std_away = 'na'

            # Additional columns used for highlighting
            if novelty == 'strain':
                tree_keep_highest_agree_higher = RANKS.index(keep_highest_agree) < RANKS.index('species')
                tree_chim_highest_agree_higher = RANKS.index(chim_highest_agree) < RANKS.index('species')
            else:
                tree_keep_highest_agree_higher = RANKS.index(keep_highest_agree) < RANKS.index(novelty)
                tree_chim_highest_agree_higher = RANKS.index(chim_highest_agree) < RANKS.index(novelty)

            # This could be a case of taxonomic inflation, need to do further checks
            # to see if this would actually inflate taxonomy
            if keep_is_chim:
                is_inflation = check_if_taxonomic_inflation(
                    gid,
                    novelty,
                    lst_tax_keep_half,
                    d_gid_to_tax,
                    d_tax_to_gids,
                    tree_keep_highest_agree_higher
                )

                tree_taxon_inflation = 'yes' if is_inflation else 'no'

            else:
                tree_taxon_inflation = 'no'

            """
            Calculate the consensus column to find out if this is truly a chimera.
            
            In this case, ANI is the source of truth, when it's not inconclusive.
            If ANI is not available, defer to the marker result.
            
            We will also assume that anything below a Z-score of 1.65 is fine.
            """
            # if n_std_away == 'na' or float(n_std_away) >= 1.65:
            if ani_keep_half_is_wrong != 'inconclusive':

                # Since the keep half is correct, we only need to do the taxon inflation check here
                if ani_keep_half_is_wrong == 'no':
                    consensus_taxon_inflation = 'no'

                # Otherwise, we can use the ANI data to determine if this is taxon inflation
                else:
                    if novelty == 'strain':
                        consensus_taxon_inflation = 'no'
                    else:
                        if novelty[0] in {x[0] for x in ani_keep_tax_lst}:
                            consensus_taxon_inflation = 'yes'
                        else:
                            consensus_taxon_inflation = 'no'

                # Best case, we can just use the ANI data
                if ani_chim_half_is_wrong != 'inconclusive':
                    consensus_chimera = 'yes' if ani_keep_half_is_wrong == 'yes' or ani_chim_half_is_wrong == 'yes' else 'no'

                # Since the chimeric half data is not available, use the tree
                else:
                    consensus_chimera = 'yes' if ani_keep_half_is_wrong == 'yes' or chim_is_chim else 'no'

            else:

                # Since taxonomic inflation relies on the keep half, we will use the tree regardless
                consensus_taxon_inflation = tree_taxon_inflation

                # In this case, we are able to use the chimeric information from FastANI (on the chimeric half).
                if ani_chim_half_is_wrong != 'inconclusive':
                    consensus_chimera = 'yes' if keep_is_chim or ani_chim_half_is_wrong == 'yes' else 'no'

                # We must rely solely on the marker information
                else:
                    consensus_chimera = 'yes' if keep_is_chim or chim_is_chim else 'no'

            # # This falls within the expected pd range
            # else:
            #     consensus_chimera = 'no'
            #     consensus_taxon_inflation = 'no'

            # Perform an override of the FastANI row data, and use the top hit information
            # this is because the current dataframe only contains the values for those with
            # sufficient alignment fraction.
            # It would be more useful to display the top hit instead of 0.
            if fastani_row['keep_ani'] == 0:
                if gid in d_gid_to_ani_raw_top_hit:
                    fastani_row['keep_ani'] = d_gid_to_ani_raw_top_hit[gid]['ani']
                    fastani_row['keep_af'] = round(d_gid_to_ani_raw_top_hit[gid]['af'], 2)
                    fastani_row['keep_sp_rep'] = d_gid_to_ani_raw_top_hit[gid]['ref']
            if fastani_row['disc_ani'] == 0:
                if f'{gid}_C' in d_gid_to_ani_raw_top_hit:
                    fastani_row['disc_ani'] = d_gid_to_ani_raw_top_hit[f'{gid}_C']['ani']
                    fastani_row['disc_af'] = round(d_gid_to_ani_raw_top_hit[f'{gid}_C']['af'], 2)
                    fastani_row['disc_sp_rep'] = d_gid_to_ani_raw_top_hit[f'{gid}_C']['ref']

            # Add the row
            rows.append({
                'gid': gid,
                'contig_count': d_gid_to_contig_count[gid],
                'taxonomy_r207': d_gid_to_tax[gid],
                'gunc_tax_level': gunc_row['taxonomic_level'],
                'gunc_css': gunc_row['clade_separation_score'],
                'gunc_rss': gunc_row['reference_representation_score'],
                'gunc_contam_prop': gunc_row['contamination_portion'],
                'gunc_db': gunc_row['source'],
                'gtdb_belongs_to_sp_where_all_fail_n_in_species': d_gid_to_sp_where_all_fail.get(gid, 'no'),
                'tree_pd_between_halves': d_gid_to_pd[gid],
                'tree_pd_z_score': n_std_away,
                'tree_nodes_between_halves': n_jumps,
                'tree_novelty': novelty,
                'tree_keep_highest_agree': keep_highest_agree,
                'tree_keep_highest_agree_HIGHER': tree_keep_highest_agree_higher,
                'tree_chim_highest_agree': chim_highest_agree,
                'tree_chim_highest_agree_HIGHER': tree_chim_highest_agree_higher,
                'tree_taxonomy_keep_half': ';'.join(lst_tax_keep_half),
                'tree_taxonomy_chim_half': ';'.join(lst_tax_chim_half),
                'tree_chimera': 'yes' if keep_is_chim or chim_is_chim else 'no',
                'tree_taxon_inflation': tree_taxon_inflation,
                'ani_identity_keep_half': fastani_row['keep_ani'],
                'ani_align_frac_keep_half': fastani_row['keep_af'],
                'ani_identity_chim_half': fastani_row['disc_ani'],
                'ani_align_frac_chim_half': fastani_row['disc_af'],
                'ani_keep_half_rep': fastani_row['keep_sp_rep'],
                'ani_chim_half_rep': fastani_row['disc_sp_rep'],
                'ani_taxonomy_keep_half': ';'.join(ani_keep_tax_lst),
                'ani_taxonomy_keep_half_full': d_gid_to_tax.get(fastani_row.get('keep_sp_rep', 'na'), 'na'),
                'ani_taxonomy_chim_half': ';'.join(ani_chim_tax_lst),
                'ani_taxonomy_chim_half_full':  d_gid_to_tax.get(fastani_row.get('disc_sp_rep', 'na'), 'na'),
                'ani_keep_half_is_wrong': ani_keep_half_is_wrong,
                'ani_chim_half_is_wrong': ani_chim_half_is_wrong,
                'consensus_chimera': consensus_chimera,
                'consensus_taxon_inflation': consensus_taxon_inflation,
            })

        df = pd.DataFrame(rows)

        df.to_csv('/tmp/test_out.tsv', sep='\t', index=False)
        return


def parse_phylorank_tax_string(taxonomy):
    out = defaultdict(list)
    for rank in taxonomy.split('; '):
        out[rank[0]].append(rank)
    return out


def calculate_highest_agree_rank(fail_gids, d_gid_to_tax_novelty, d_gid_to_tax, d_gid_to_tree_tax):
    out = dict()

    for gid in tqdm(sorted(fail_gids)):

        gid_true_tax = d_gid_to_tax[gid]
        gid_true_tax_split = gid_true_tax.split(';')
        gid_novelty = d_gid_to_tax_novelty[gid]
        gid_novelty_idx = D_TAX_NOVELTY_TO_INDEX[gid_novelty]

        for suffix in ('', '_C'):
            gid_suffix = f'{gid}{suffix}'
            # gid_split_tax = d_gid_to_split_tax[f'TEST_{gid_suffix}']
            # gid_split_tax_parsed = parse_phylorank_tax_string(gid_split_tax)
            gid_split_tax_parsed = d_gid_to_tree_tax[gid_suffix]

            """Three main cases:
            novel at class level
            true: d__A; p__B; c__C; o__D; f__E; g__F; s__G
            test: d__A; p__B; c__;  o__;   f__;   g__;   s__
            If the test taxonomy is blank, then need to go up to find a rank that agrees.
            
            true: d__A; p__B; c__C; o__D; f__E; g__F; s__G
            test: d__A; p__B; c__D;  o__;   f__;   g__;   s__
            If the test taxonomy is disagrees, then go up.
            
          true: d__A; p__B; c__C; o__D; f__E; g__F; s__G
            test: d__A; p__B; c__;  o__;   f__;   g__;   s__
            If the test taxonomy agrees, then go down
            """
            # Do a quick pass over all ranks and do a check if they agree
            idx_to_rank_correct = list()
            for rank_idx, true_rank in enumerate(gid_true_tax_split):
                lst_test_rank = gid_split_tax_parsed[true_rank[0]]
                if len(lst_test_rank) == 0:
                    idx_to_rank_correct.append(False)
                elif len(lst_test_rank) != 1:
                    # If its the first rank, then it's fine (note the order is reversed due to how the tree tax is parsed)
                    idx_to_rank_correct.append(true_rank == lst_test_rank[0])
                    # idx_to_rank_correct.append(False)
                else:
                    if rank_idx >= gid_novelty_idx:
                        idx_to_rank_correct.append(
                            true_rank == lst_test_rank[0] or lst_test_rank[0] == f'{true_rank[0]}__')
                    elif rank_idx == gid_novelty_idx - 1:
                        idx_to_rank_correct.append(
                            true_rank == lst_test_rank[0] and lst_test_rank[0] != f'{true_rank[0]}__')
                    else:
                        idx_to_rank_correct.append(true_rank == lst_test_rank[0])

            # Do a quick check to see if we can skip processing this genome
            if all(idx_to_rank_correct):
                out[gid_suffix] = 'species'

            # If the genome is still correct at the correct taxonomic novelty rank, then
            # iterate downwards to check if the other ranks are correct
            elif idx_to_rank_correct[gid_novelty_idx - 1]:
                best_rank_idx = gid_novelty_idx - 1
                for rank_idx in range(gid_novelty_idx, len(gid_true_tax_split)):
                    is_correct = idx_to_rank_correct[rank_idx]
                    if not is_correct:
                        break
                    best_rank_idx = rank_idx
                out[gid_suffix] = RANKS[best_rank_idx]

            # Otherwise, we need to go up and find which rank is first congruent
            else:
                best_rank_idx = -1
                for rank_idx in reversed(range(gid_novelty_idx)):
                    is_correct = idx_to_rank_correct[rank_idx]
                    if is_correct:
                        best_rank_idx = rank_idx
                        break
                out[gid_suffix] = RANKS[best_rank_idx]

    return out


def get_taxonomy_from_tree(tree):
    out = dict()

    for leaf_node in tqdm(tree.leaf_node_iter(), total=len(tree.taxon_namespace)):

        # Only interested in those that are specified
        if not leaf_node.taxon.label.startswith('TEST_'):
            continue

        cur_tax = defaultdict(list)

        cur_node = leaf_node

        while cur_node is not None:
            _, taxon, _ = parse_label(cur_node.label)

            if taxon is not None:
                for rank in taxon.split('; '):
                    cur_tax[rank[0]].append(rank)
            cur_node = cur_node.parent_node

        # Append any missing ranks
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            if rank not in cur_tax:
                cur_tax[rank].append(f'{rank}__')

        # Save the taxonomy to the output dictionary
        out[leaf_node.taxon.label.replace('TEST_', '')] = cur_tax
    return out



def calculate_for_strain(df_out, novelty, min_value=None, max_value=None, plot_min=None, plot_max=None, plot=True):
    df_strain = df_out[df_out['novelty'] == novelty]
    if min_value:
        df_strain = df_strain[df_strain['pd'] >= min_value]
    if max_value:
        df_strain = df_strain[df_strain['pd'] <= max_value]

    if plot:
        sns.histplot(df_strain, x='pd', log_scale=True, hue='replicate_id')
        plt.show()


    strain_values = list()
    for i in range(N_MARKER_SPLIT_TRUE):

        if plot:
            fig, ax = plt.subplots(1, 1)

        r = df_strain[df_strain['replicate_id'] == i]['pd'].values
        r = np.log(r)
        loc, scale = norm.fit(r)

        median = norm.median(loc, scale)
        mean = norm.mean(loc, scale)
        std = lognorm.std(loc, scale)

        strain_values.append({
            'replicate_id': i,
            'loc': loc,
            'scale': scale,
            'median': median,
            'mean': mean,
            'std': std,
        })

        if plot:
            x = np.linspace(0, r.max() * 2, 10000)
            ax.plot(x, norm.pdf(x, loc, scale), 'r-', lw=1, alpha=0.6, label=f'norm pdf rep {i}')
            sns.kdeplot(r, ax=ax)
            ax.legend(loc='best', frameon=False)
            plt.show()

    df_strain_values = pd.DataFrame(strain_values)
    strain_loc = df_strain_values['loc'].mean()
    strain_scale = df_strain_values['scale'].mean()

    if plot:
        fig, ax = plt.subplots(1, 1)
        r = df_strain['pd'].values

        x = np.linspace(0, 1000, 10000)
        ax.plot(x, norm.pdf(x, strain_loc, strain_scale), 'r-', lw=1, alpha=0.6)
        for i in range(N_MARKER_SPLIT_TRUE):
            r = np.log(df_strain[df_strain['replicate_id'] == i]['pd'].values)
            sns.kdeplot(r, ax=ax, alpha=0.2)
        # sns.kdeplot(np.log(df_strain['pd'].values), ax=ax)
        if plot_min is not None and plot_max is not None:
            ax.set_xlim(plot_min, plot_max)
        plt.show()

    avg_mean = norm.mean(strain_loc, strain_scale)
    avg_std = norm.std(strain_loc, strain_scale)

    print(f'Novelty: {novelty}, loc={strain_loc}, scale={strain_scale}, mean={avg_mean:.4f}, std={avg_std:.4f}')
    return float(strain_loc), float(strain_scale), float(avg_mean), float(avg_std)

def calculate_for_species(df_out):
    cur_rank = 'species'
    df_rank = df_out[(df_out['novelty'] == cur_rank)]
    sns.kdeplot(df_rank, x='pd', hue='replicate_id')
    plt.show()


    out_values = list()
    for i in range(N_MARKER_SPLIT_TRUE):
        fig, ax = plt.subplots(1, 1)

        r = df_rank[df_rank['replicate_id'] == i]['pd'].values
        s, loc, scale = lognorm.fit(r, scale=100, loc=-5)

        median = lognorm.median(s, loc, scale)
        mean = lognorm.mean(s, loc, scale)
        std = lognorm.std(s, loc, scale)

        out_values.append({
            'replicate_id': i,
            's': s,
            'loc': loc,
            'scale': scale,
            'median': median,
            'mean': mean,
            'std': std,
        })

        x = np.linspace(r.min(), r.max(), 10000)
        ax.plot(x, lognorm.pdf(x, s, loc, scale), 'r-', lw=1, alpha=0.6, label=f'lognorm pdf {i}')
        sns.kdeplot(df_rank[df_rank['replicate_id'] == i], x='pd', ax=ax)
        ax.legend(loc='best', frameon=False)
        plt.show()

    df_strain_values = pd.DataFrame(out_values)
    strain_s = df_strain_values['s'].mean()
    strain_loc = df_strain_values['loc'].mean()
    strain_scale = df_strain_values['scale'].mean()

    fig, ax = plt.subplots(1, 1)
    r = df_rank['pd'].values

    x = np.linspace(r.min(), r.max(), 10000)
    ax.plot(x, lognorm.pdf(x, strain_s, strain_loc, strain_scale), 'r-', lw=1, alpha=0.6)
    sns.kdeplot(df_rank, x='pd', ax=ax)
    ax.legend(loc='best', frameon=False)
    plt.show()
    print()

    return strain_s, strain_loc, strain_scale


def calculate_for_genus(df_out):
    cur_rank = 'genus'
    df_rank = df_out[(df_out['novelty'] == cur_rank)]
    sns.kdeplot(df_rank, x='pd', hue='replicate_id')
    plt.show()


    out_values = list()
    for i in range(N_MARKER_SPLIT_TRUE):
        fig, ax = plt.subplots(1, 1)

        r = df_rank[df_rank['replicate_id'] == i]['pd'].values
        s, loc, scale = lognorm.fit(r, scale=100, loc=-5)

        median = lognorm.median(s, loc, scale)
        mean = lognorm.mean(s, loc, scale)
        std = lognorm.std(s, loc, scale)

        out_values.append({
            'replicate_id': i,
            's': s,
            'loc': loc,
            'scale': scale,
            'median': median,
            'mean': mean,
            'std': std,
        })

        x = np.linspace(r.min(), r.max(), 10000)
        ax.plot(x, lognorm.pdf(x, s, loc, scale), 'r-', lw=1, alpha=0.6, label=f'lognorm pdf {i}')
        sns.kdeplot(df_rank[df_rank['replicate_id'] == i], x='pd', ax=ax)
        ax.legend(loc='best', frameon=False)
        plt.show()

    df_strain_values = pd.DataFrame(out_values)
    strain_s = df_strain_values['s'].mean()
    strain_loc = df_strain_values['loc'].mean()
    strain_scale = df_strain_values['scale'].mean()

    fig, ax = plt.subplots(1, 1)
    r = df_rank['pd'].values

    x = np.linspace(r.min(), r.max(), 10000)
    ax.plot(x, lognorm.pdf(x, strain_s, strain_loc, strain_scale), 'r-', lw=1, alpha=0.6)
    sns.kdeplot(df_rank, x='pd', ax=ax)
    ax.legend(loc='best', frameon=False)
    plt.show()
    print()

    return strain_s, strain_loc, strain_scale


def calculate_for_family(df_out):
    # There are too few points for the family level.


    cur_rank = 'family'
    df_rank = df_out[(df_out['novelty'] == cur_rank)]
    sns.kdeplot(df_rank, x='pd', hue='replicate_id')
    plt.show()


    out_values = list()
    for i in range(N_MARKER_SPLIT_TRUE):
        fig, ax = plt.subplots(1, 1)

        r = df_rank[df_rank['replicate_id'] == i]['pd'].values
        loc, scale = norm.fit(r, scale=100, loc=-5)

        median = norm.median(loc, scale)
        mean = norm.mean(loc, scale)
        std = norm.std(loc, scale)

        out_values.append({
            'replicate_id': i,
            'loc': loc,
            'scale': scale,
            'median': median,
            'mean': mean,
            'std': std,
        })

        x = np.linspace(r.min(), r.max(), 10000)
        ax.plot(x, norm.pdf(x, loc, scale), 'r-', lw=1, alpha=0.6, label=f'norm pdf {i}')
        sns.kdeplot(df_rank[df_rank['replicate_id'] == i], x='pd', ax=ax)
        ax.legend(loc='best', frameon=False)
        plt.show()

    df_strain_values = pd.DataFrame(out_values)
    strain_s = df_strain_values['s'].mean()
    strain_loc = df_strain_values['loc'].mean()
    strain_scale = df_strain_values['scale'].mean()

    fig, ax = plt.subplots(1, 1)
    r = df_rank['pd'].values

    x = np.linspace(r.min(), r.max(), 10000)
    ax.plot(x, norm.pdf(x, strain_s, strain_loc, strain_scale), 'r-', lw=1, alpha=0.6)
    sns.kdeplot(df_rank, x='pd', ax=ax)
    ax.legend(loc='best', frameon=False)
    plt.show()
    print()

    return strain_s, strain_loc, strain_scale


def check_if_taxonomic_inflation(gid, novelty, lst_tax_keep, d_gid_to_tax, d_tax_to_gid, tree_keep_highest_agree_higher):

    """
    The most trivial case is when a genome is the only representative (e.g. novel at the class level)
    and the resulting taxonomy puts this genome in a new class.
    """
    if tree_keep_highest_agree_higher:

        if novelty == 'strain':
            """
              If this genome is novel at the strain level, then by definition
              it will not cause inflation by leaving the original group.
              However, if it causes a novel group it will cause inflation.
              """
            return any(len(x) == 3 for x in lst_tax_keep)

        else:
            """Otherwise, since this genome disagrees at a higher rank than the
            taxonomic novelty level. Then this would have left the original group,
            i.e. a loss of that taxonomic rank. This would be deflation."""
            return True

    raise Exception('??')

def get_species_where_all_fail(df_css, df_meta):
    out = dict()

    fail_gids = set(df_css.index)

    d_species_to_gids = defaultdict(set)
    for row in df_meta.itertuples():
        d_species_to_gids[row.species].add(row.Index)

    for species, gids in d_species_to_gids.items():
        if gids.issubset(fail_gids):
            for gid in gids:
                out[gid] = len(gids)
    return out
