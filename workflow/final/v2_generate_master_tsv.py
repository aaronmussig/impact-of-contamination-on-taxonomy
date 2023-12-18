import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import lognorm
from scipy.stats import norm
from tqdm import tqdm

from workflow.config import DIR_OUT_FINAL
from workflow.external.gtdb_metadata import GtdbMetadataR207, GtdbMetadataR207Full
from workflow.fastani_contig_split.b_run_fastani import FastAniContigSplitRunFastAni
from workflow.fastani_contig_split.d_report_results_relaxed import FastAniContigSplitReportResultsFastAniRelaxed
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetTsv
from workflow.util.log import log
from workflow.util.taxonomy import calculate_taxonomic_novelty, D_TAX_NOVELTY_TO_INDEX
from workflow.util.tree import parse_label
from workflow.v2_fasttree_marker_split.e_analyse_decorated import V2FastTreeMarkerSplitAnalyseDecoratedFastTree
from workflow.v2_pplacer_marker_split.c_get_taxonomy_from_tree import V2PplacerMarkerSplitGetTaxonomyFromTree

RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')

RANKS_STRAIN = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')

N_MARKER_SPLIT_TRUE = 10

D_RANK_PREFIX_TO_RANK = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}

RANK_PREFIX_TO_IDX = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}

class V2FinalGenerateMasterTsv(LuigiTask):

    def requires(self):
        req = {
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
            'meta_full': GtdbMetadataR207Full(),
            'fastani_split': FastAniContigSplitReportResultsFastAniRelaxed(),
            'ani_raw': FastAniContigSplitRunFastAni(),
            'v2_marker_split_rep_tree': V2FastTreeMarkerSplitAnalyseDecoratedFastTree(remove_pct=50),
            'pplacer': V2PplacerMarkerSplitGetTaxonomyFromTree(remove_pct=50, n_trees_output=8),
        }
        return req

    def output(self):
        return LocalTargetTsv(os.path.join(DIR_OUT_FINAL, 'master.tsv'))

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

    def load_type_material(self):
        out = dict()
        df = self.input()['meta_full'].maybe_read_cached()
        for gid, cat in df['ncbi_genome_category'].to_dict().items():
            if cat == 'derived from metagenome':
                out[gid] = 'MAG'
            elif cat == 'derived from single cell':
                out[gid] = 'SAG'
            else:
                out[gid] = 'isolate'
        return out

    def run(self):
        log(f'V2FinalGenerateMasterTsv', title=True)
        self.make_output_dirs()

        log('Loading pplacer taxonomy (non-reps)')
        df_pplacer = self.input()['pplacer'].maybe_read_cached()
        df_pplacer.set_index('gid', inplace=True)
        d_pplacer_gid_to_tax = df_pplacer['tax'].to_dict()

        log('Loading analysis of marker rep tree split')
        df_marker_split_rep_tree = self.input()['v2_marker_split_rep_tree'].maybe_read_cached()
        df_marker_split_rep_tree.set_index('gid', inplace=True)

        log('Loading fastani split results')
        df_fastani = self.input()['fastani_split'].maybe_read_cached()
        df_fastani.set_index('gid', inplace=True)

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()
        d_gid_to_contig_count = df_meta['contig_count'].to_dict()
        s_bac_gids = frozenset(df_meta[df_meta['domain'] == 'd__Bacteria'].index)
        s_bac_rep_gids = frozenset(
            df_meta[(df_meta['domain'] == 'd__Bacteria') & (df_meta['gtdb_representative'] == 't')].index)

        d_sp_to_rep_gid = dict()
        for gid in s_bac_rep_gids:
            d_sp_to_rep_gid[d_gid_to_tax[gid].split(';')[-1]] = gid

        d_sp_to_gids = defaultdict(set)
        for gid, tax in d_gid_to_tax.items():
            d_sp_to_gids[tax.split(';')[-1]].add(gid)

        log('Loading MaxCSS')
        df_css = self.input()['max_css'].maybe_read_cached()
        fail_gids = frozenset(df_css.index).intersection(s_bac_gids)

        log('Getting taxonomic novelty')
        d_gid_to_tax_novelty, d_tax_novelty_to_gid = calculate_taxonomic_novelty(d_gid_to_tax)
        log('Getting count of failed genomes in each species')
        d_gid_to_n_failed_gids_in_sp, s_gids_where_in_sp_with_all_fail = get_count_of_failed_gids_in_sp(df_css, df_meta)

        log('Loading genome source material')
        d_gid_to_source = self.load_type_material()

        log('Collecting raw FastANI data')
        d_gid_to_ani_raw_top_hit = self.load_top_hit_from_ani_raw_data()

        log('Creating dataframe rows')
        rows = list()
        for gid in tqdm(sorted(fail_gids)):

            gunc_row = df_css.loc[gid]
            novelty = d_gid_to_tax_novelty[gid]
            fastani_row = df_fastani.loc[gid]
            cur_sp = d_gid_to_tax[gid].split(';')[-1]
            cur_sp_rep_gid = d_sp_to_rep_gid[cur_sp]
            gid_tax = d_gid_to_tax[gid]
            gid_tax_split = gid_tax.split(';')
            gid_c = f'{gid}_C'

            base_row = {
                'gid': gid,
                'contig_count': d_gid_to_contig_count[gid],
                'gtdb_source': d_gid_to_source[gid],
                'gtdb_taxonomy': gid_tax,
                'gtdb_is_sp_rep': 'yes' if gid in s_bac_rep_gids else 'no',
                'gtdb_novelty': novelty,
                'gtdb_n_genomes_in_sp_cluster': len(d_sp_to_gids[cur_sp]),
                'n_gids_in_sp_fail': d_gid_to_n_failed_gids_in_sp[gid],
                'all_gids_in_sp_fail': 'yes' if gid in s_gids_where_in_sp_with_all_fail else 'no',
                'sp_rep_pass_or_fail': 'fail' if cur_sp_rep_gid in fail_gids else 'pass',
                'gunc_tax_level': gunc_row['taxonomic_level'],
                'gunc_css': gunc_row['clade_separation_score'],
                'gunc_rss': gunc_row['reference_representation_score'],
                'gunc_contam_prop': gunc_row['contamination_portion'],
                'gunc_db': gunc_row['source'],
                'tree_rank_highest_agree': 'todo',  # done
                'tree_taxonomy_disagree': 'todo',  # done
                'tree_congruent': 'todo',  # done
                'tree_taxon_inflation': 'todo',
                'ani_rank_highest_agree': 'todo',  # done
                'ani_identity': 'todo',  # done
                'ani_af': 'todo',  # done
                'ani_taxonomy_disagree': 'todo',  # done
                'ani_representative': 'todo',  # done
                'ani_congruent': 'todo',  # done
                'ani_taxon_inflation': 'todo',
                'ani_chim_rank_highest_agree': 'todo',
                'ani_chim_identity': 'todo',
                'ani_chim_af': 'todo',
                'ani_chim_taxonomy_disagree': 'todo',
                'ani_chim_representative': fastani_row['disc_sp_rep'],
                'ani_chim_congruent': 'todo',
                'consensus_congruent': 'todo',
                'consensus_taxon_inflation': 'todo',
            }

            # Append the chimeric information (mainly for debugging information)
            if gid_c in d_gid_to_ani_raw_top_hit:
                chim_tax = d_gid_to_tax[d_gid_to_ani_raw_top_hit[gid_c]['ref']]
                chim_ani = d_gid_to_ani_raw_top_hit[gid_c]['ani']
                chim_af = round(d_gid_to_ani_raw_top_hit[gid_c]['af'], 4)
                base_row['ani_chim_rank_highest_agree'] = get_highest_agreed_rank(gid_tax, chim_tax)
                base_row['ani_chim_identity'] = chim_ani
                base_row['ani_chim_af'] = chim_af
                base_row['ani_chim_taxonomy_disagree'] = trim_ranks_that_agree(gid_tax, chim_tax)
                base_row['ani_chim_representative'] = d_gid_to_ani_raw_top_hit[gid_c]['ref']
                base_row['ani_chim_congruent'] = 'yes' if base_row['ani_chim_rank_highest_agree'] == 'species' else 'no'

                if chim_ani >= 95 and chim_af >= 0.5:
                    if gid_tax == chim_tax:
                        base_row['ani_chim_congruent'] = 'yes'
                    else:
                        base_row['ani_chim_congruent'] = 'no'
                else:
                    base_row['ani_chim_congruent'] = 'inconclusive'

            else:
                base_row['ani_chim_rank_highest_agree'] = ''
                base_row['ani_chim_identity'] = ''
                base_row['ani_chim_af'] = ''
                base_row['ani_chim_taxonomy_disagree'] = ''
                base_row['ani_chim_representative'] = ''
                base_row['ani_chim_congruent'] = 'inconclusive'

            # If we have a hit <95% ANI use the raw data
            if fastani_row['keep_ani'] == 0:
                base_row['ani_congruent'] = 'inconclusive'
                base_row['ani_rank_highest_agree'] = ''
                base_row['ani_taxon_inflation'] = ''

                # We can use the top hit data that has <95% ani / <0.5 AF
                if gid in d_gid_to_ani_raw_top_hit:
                    base_row['ani_identity'] = d_gid_to_ani_raw_top_hit[gid]['ani']
                    base_row['ani_af'] = round(d_gid_to_ani_raw_top_hit[gid]['af'], 4)
                    base_row['ani_taxonomy_disagree'] = trim_ranks_that_agree(
                        gid_tax,
                        d_gid_to_tax[d_gid_to_ani_raw_top_hit[gid]['ref']]
                    )
                    base_row['ani_representative'] = d_gid_to_ani_raw_top_hit[gid]['ref']
                else:
                    base_row['ani_identity'] = ''
                    base_row['ani_af'] = ''
                    base_row['ani_taxonomy_disagree'] = ''
                    base_row['ani_representative'] = ''

            # We have a hit within >95% ANI, use that data instead
            else:
                base_row['ani_identity'] = fastani_row['keep_ani']
                base_row['ani_af'] = fastani_row['keep_af']
                base_row['ani_taxonomy_disagree'] = trim_ranks_that_agree(gid_tax, fastani_row['keep_tax'])
                base_row['ani_representative'] = fastani_row['keep_sp_rep']
                base_row['ani_rank_highest_agree'] = get_highest_agreed_rank(gid_tax, fastani_row['keep_tax'])

                # rep gids can't compare to themselves
                # if gid in s_bac_rep_gids:
                #     base_row['ani_congruent'] = 'inconclusive'
                #     base_row['ani_taxon_inflation'] = ''
                # else:
                if fastani_row['keep_ani'] >= 95 and fastani_row['keep_af'] >= 0.5:
                    if gid_tax == fastani_row['keep_tax']:
                        base_row['ani_congruent'] = 'yes'
                    else:
                        base_row['ani_congruent'] = 'no'
                else:
                    base_row['ani_congruent'] = 'inconclusive'

                if novelty != 'strain':
                    base_row['ani_taxon_inflation'] = 'todo'
                else:
                    base_row['ani_taxon_inflation'] = 'no'

            # append the tree data
            if gid in s_bac_rep_gids:
                tree_row = df_marker_split_rep_tree.loc[f'TEST_{gid}']

                if tree_row['agree'] == True:
                    base_row['tree_rank_highest_agree'] = ''
                    base_row['tree_taxonomy_disagree'] = ''
                    base_row['tree_congruent'] = 'yes'
                    base_row['tree_taxon_inflation'] = 'no'
                else:
                    base_row['tree_rank_highest_agree'] = D_RANK_PREFIX_TO_RANK[tree_row['highest_disagree_rank']]
                    base_row['tree_taxonomy_disagree'] = trim_ranks_that_agree(gid_tax, tree_row['test_tax'])
                    base_row['tree_congruent'] = 'no'

                    contains_novel_ranks = any([len(x) == 3 for x in tree_row['test_tax'].split(';')])
                    if novelty == 'strain':
                        base_row['tree_taxon_inflation'] = 'yes' if contains_novel_ranks else 'no'
                    else:
                        base_row['tree_taxon_inflation'] = 'yes'

                    # Taxon inflation for the tree can either be: the half has jumped to another taxon
                    # this is inflation (if it is not novel at the strain level)

                    # or, if it's creating a new taxon (e.g. s__)

            # this was a non-rep gid so take the results from pplacer
            # note, as this is a non-rep gid it will be a strain level novelty by definition
            else:
                pplacer_tax = d_pplacer_gid_to_tax[gid]
                pplacer_tax_split = pplacer_tax.split(';')
                ref_tax_without_sp = ';'.join(gid_tax_split[0:6])
                pplacer_tax_without_sp = ';'.join(pplacer_tax_split[0:6])
                is_tree_congruent = ref_tax_without_sp == pplacer_tax_without_sp

                if is_tree_congruent:
                    base_row['tree_rank_highest_agree'] = ''
                    base_row['tree_taxonomy_disagree'] = ''
                    base_row['tree_congruent'] = 'yes'
                    base_row['tree_taxon_inflation'] = 'no'
                else:
                    base_row['tree_rank_highest_agree'] = get_highest_agreed_rank(gid_tax, pplacer_tax)
                    base_row['tree_taxonomy_disagree'] = trim_ranks_that_agree(gid_tax, pplacer_tax)
                    base_row['tree_congruent'] = 'no'
                    base_row['tree_taxon_inflation'] = 'todo'
                    raise Exception('TODO')




            # Save the row
            rows.append(base_row)

        df = pd.DataFrame(rows)

        df.to_csv('/tmp/test_out_v2.tsv', sep='\t', index=False)
        return


def parse_phylorank_tax_string(taxonomy):
    out = defaultdict(list)
    for rank in taxonomy.split(';'):
        out[rank[0]].append(rank)
    return out

def trim_ranks_that_agree(true_tax, test_tax):
    d_tax_true = {x[0]: x for x in true_tax.split(';')}

    d_test_tax = parse_phylorank_tax_string(test_tax)
    out = list()
    append_remaining = False
    for i, rank in enumerate('dpcofgs'):
        true_rank = d_tax_true[rank]
        if len(d_test_tax[rank]) == 1:
            test_rank = d_test_tax[rank][0]

        elif len(d_test_tax[rank]) > 1:

            # if it's a nested rank then it's ok
            if true_rank in d_test_tax[rank]:
                test_rank = true_rank
            else:
                test_rank = d_test_tax[rank][0]

        else:
            raise Exception('?')

        if true_rank != test_rank:
            append_remaining = True

        if append_remaining:
            out.extend(d_test_tax[rank])

    if len(out) == 0:
        return ''
    else:
        return ';'.join(out)


def get_highest_agreed_rank(true_tax, test_tax):
    d_tax_true = {x[0]: x for x in true_tax.split(';')}

    d_test_tax = parse_phylorank_tax_string(test_tax)
    all_ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    for i, rank in enumerate(all_ranks):
        true_rank = d_tax_true[rank]
        if len(d_test_tax[rank]) == 1:
            test_rank = d_test_tax[rank][0]

        elif len(d_test_tax[rank]) > 1:

            # if it's a nested rank then it's ok
            if true_rank in d_test_tax[rank]:
                test_rank = true_rank
            else:
                test_rank = d_test_tax[rank][0]

        else:
            raise Exception('?')

        if true_rank != test_rank:
            return D_RANK_PREFIX_TO_RANK[all_ranks[i-1]]

    return 'species'


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


def check_if_taxonomic_inflation(gid, novelty, lst_tax_keep, d_gid_to_tax, d_tax_to_gid,
                                 tree_keep_highest_agree_higher):
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


def get_count_of_failed_gids_in_sp(df_css, df_meta):
    out = dict()

    fail_gids = set(df_css.index)

    d_species_to_gids = defaultdict(set)
    for row in df_meta.itertuples():
        d_species_to_gids[row.species].add(row.Index)

    out_all_fail = set()

    for species, gids in d_species_to_gids.items():
        for gid in gids:
            out[gid] = len(gids.intersection(fail_gids))
        if gids.issubset(fail_gids):
            for gid in gids:
                out_all_fail.add(gid)
    return out, out_all_fail
