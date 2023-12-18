import multiprocessing as mp
import os
from collections import Counter, defaultdict
from typing import Set, Dict

import dendropy
import luigi
import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, DEBUG, DIR_OUT_FASTTREE_MARKER_SPLIT
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_tree import GtdbTreeBacR207
from workflow.fasttree_full_tree_and_failed_no_trim.e_decorate_fasttree import FastTreeFullTreeAndFailedNoTrimDecorate
from workflow.fasttree_full_tree_non_reps.gunc_c_decorate_fasttree import \
    FastTreeFullSetNonRepsCreateBatchesGuncDecorate
from workflow.fasttree_marker_analysis.random_d_analyse_decorated import decorated_tax_to_agreed, get_taxon_to_gids
from workflow.fasttree_marker_split.e_decorate_fasttree import FastTreeMarkerSplitDecorateFastTree
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from workflow.util.tree import parse_label, get_tree_node_to_desc_taxa


def read_tax_file(path):
    d_gid_to_decorated_tax = dict()
    with open(path) as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            d_gid_to_decorated_tax[gid] = tax.replace('; ', ';')
    return d_gid_to_decorated_tax


def read_gids_file(path):
    out = set()
    with open(path) as f:
        for line in f.readlines():
            if line.strip() == '':
                continue
            out.add(line.strip())
    return out


class FastTreeMarkerSplitDecorateAnalyseFastTree(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_MARKER_SPLIT, f'pct_{self.target_pct}')

    def requires(self):
        return {
            'decorated': FastTreeMarkerSplitDecorateFastTree(target_pct=self.target_pct),
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'r207_bac_ref_tree': GtdbTreeBacR207(),
            'no_trim_tree': FastTreeFullTreeAndFailedNoTrimDecorate(),
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, 'decorated_analysis_results.tsv'))

    def run(self):
        log(f'FastTreeMarkerSplitDecorateAnalyseFastTree(pct={self.target_pct})', title=True)
        self.make_output_dirs()

        # Set paths for the trimmed tree
        trim_root_dir = os.path.dirname(self.input()['decorated'].path)
        trim_decorated_tax = os.path.join(trim_root_dir, 'decorated.tree-taxonomy')
        trim_d_gid_to_tax = read_tax_file(trim_decorated_tax)

        # Set the paths for the untrimmed tree
        untrim_root_dir = os.path.dirname(self.input()['no_trim_tree'].path)
        untrim_decorated_tax = os.path.join(untrim_root_dir, 'decorated.tree-taxonomy')
        untrim_d_gid_to_tax = read_tax_file(untrim_decorated_tax)

        log('Loading untrimmed tree')
        tree_untrimmed = self.input()['no_trim_tree'].read()

        log('Loading trimmed tree')
        tree_trimmed = self.input()['decorated'].read()

        gids_where_halves_differ = get_chimeric_halves_where_tax_strings_differ(trim_d_gid_to_tax, untrim_d_gid_to_tax, tree_trimmed, tree_untrimmed)



        log('Loading untrimmed tree')
        tree_untrimmed = self.input()['no_trim_tree'].read()

        log('Loading trimmed tree')
        tree_trimmed = self.input()['decorated'].read()

        # Compare the tax strings
        compare_tax_string_dicts(trim_d_gid_to_tax, untrim_d_gid_to_tax, tree_trimmed, tree_untrimmed)
        

        """
        TODO: Get the tax string from the decorated tree, compare it with the untrimmed.
        Any that are not literal matches are chimeric?
        
        """






        # log('Loading BAC120 R207 reference tree')
        # r207_bac_ref_tree = self.input()['r207_bac_ref_tree'].read()
        # for leaf_node in r207_bac_ref_tree.leaf_node_iter():
        #     leaf_node.taxon.label = leaf_node.taxon.label[3:]

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()

        log('Loading MaxCSS file')
        df_css = self.input()['max_css'].maybe_read_cached()
        fail_gids = set(df_css.index).intersection(set(df_meta[df_meta['domain'] == 'd__Bacteria'].index))

        fail_gids_with_prefix_and_suffix = set()
        for gid in fail_gids:
            fail_gids_with_prefix_and_suffix.add(f'TEST_{gid}')
            fail_gids_with_prefix_and_suffix.add(gid)
            fail_gids_with_prefix_and_suffix.add(f'TEST_{gid}_C')



        match_nodes(tree_untrimmed, tree_trimmed, fail_gids, d_gid_to_tax)


        untrimmed_sister_taxa = get_sister_taxa(tree_untrimmed, fail_gids_with_prefix_and_suffix)
        trimmed_sister_taxa = get_sister_taxa(tree_trimmed, fail_gids_with_prefix_and_suffix)


        log('Comparing sister taxa')
        for trimmed_taxon, cur_trimmed_sister_taxa in trimmed_sister_taxa.items():
            cur_gid = trimmed_taxon.replace('TEST_', '').replace('_C', '')

            untrimmed_taxa = untrimmed_sister_taxa[cur_gid]
            trimmed_taxa = {x.replace('TEST_', '').replace('_C', '') for x in cur_trimmed_sister_taxa}

            common_gids = untrimmed_taxa.intersection(trimmed_taxa) - {cur_gid}

            print()

        log('Loading decorated tree')
        # decorated_tree = self.input()['decorated'].read()

        # Set paths
        phylo_dir = os.path.dirname(self.input()['decorated'].path)
        root_dir = os.path.dirname(phylo_dir)


        d_gid_to_taxonomy = df_meta['gtdb_taxonomy'].to_dict()

        path_decorated_tax = os.path.join(phylo_dir, 'decorated.tree-taxonomy')
        d_gid_to_decorated_tax = read_tax_file(path_decorated_tax)

        log('Getting taxonomic novelty (reps only)')
        d_taxon_to_gids = get_taxon_to_gids(df_meta[df_meta['gtdb_representative'] == 't']['gtdb_taxonomy'].to_dict())

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()

        # For each half, find the tax string of it's counterpart
        for gid in tqdm(df_max_css.index):
            if not d_gid_to_taxonomy[gid].startswith('d__Bacteria'):
                continue

            keep_tax = d_gid_to_decorated_tax[f'TEST_{gid}']
            disc_tax = d_gid_to_decorated_tax[f'TEST_{gid}_C']
            if keep_tax.split(';')[0:5] != disc_tax.split(';')[0:5]:
                print('?')

        return



        correct = set()
        wrong = set()
        for test_gid, decorated_tax in tqdm(sorted(d_gid_to_decorated_tax.items())):
            if test_gid.startswith('TEST_'):
                if test_gid.endswith('_C'):
                    gid = test_gid[5:-2]
                else:
                    gid = test_gid[5:]
            else:
                continue
            true_tax = d_gid_to_taxonomy[gid]
            true_tax_no_sp = true_tax[:true_tax.index('s__') - 1]
            if not decorated_tax.startswith(true_tax_no_sp):
            # if true_tax != decorated_tax:

                # d_decorated_tax_lst = decorated_tax_to_dict(decorated_tax)
                d_decorated_tax_suffix = decorated_tax_to_agreed(decorated_tax)

                # 1. The highest taxon (that has more than one genome) must be the same
                # 2. All descendants in that taxon, must agree on all lower ranks

                # Get the highest rank that should be taxonomically novel
                highest_rank = None
                for taxon in true_tax.split(';'):
                    taxon_count = len(d_taxon_to_gids[taxon])
                    if taxon_count == 1:
                        break
                    highest_rank = taxon

                # Must be comparing a higher rank than the singleton
                if highest_rank not in d_decorated_tax_suffix:
                    wrong.add(test_gid)
                else:
                    if not true_tax.endswith(d_decorated_tax_suffix):
                        wrong.add(test_gid)
                    else:
                        correct.add(test_gid)

                # print('mismatch', test_gid, true_tax, decorated_tax)
            else:
                correct.add(test_gid)

        for gid in sorted(wrong):
            wrong_tax = d_gid_to_decorated_tax[f'TEST_{gid}']
            true_tax = d_gid_to_taxonomy[gid]
            is_rep = df_meta.loc[gid, 'gtdb_representative'] == 't'
            expected_sp = df_meta.loc[gid, 'species']
            n_gids_in_sp = df_meta[df_meta['species'] == expected_sp].shape[0]
            print(f'{gid} rep={is_rep} ({n_gids_in_sp:,} in cluster)\n{true_tax}\n{wrong_tax}')
            print()
            print()

        # Get the equivalent node in the true reference tree
        dec_sister_taxa = get_sister_taxa(decorated_tree, {f'TEST_{x}' for x in wrong})
        ref_sister_taxa = get_sister_taxa(r207_bac_ref_tree, wrong)

        # Compare the results
        gids_to_check = set()
        for gid in wrong:
            ref_taxa = ref_sister_taxa[gid]
            dec_taxa = dec_sister_taxa[f'TEST_{gid}']
            dec_taxa_strip = {x.replace('TEST_', '') for x in dec_taxa}
            if ref_taxa != dec_taxa_strip:
                gids_to_check.add(gid)

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()

        log('Creating queue')
        queue = list()
        for gid in s_all_gids:
            row = df_max_css.loc[gid]
            queue.append((gid, row['source'], row['taxonomic_level'], self.congruence, self.target_pct))

        log('Processing queue')
        with mp.Pool(processes=int(mp.cpu_count() * 0.8)) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        out = list()
        for d_result in results:
            cur_gid = d_result['gid']
            cur_result = {
                'gid': cur_gid,
                'congruence': self.congruence,
                'target_pct': self.target_pct,
                'pct_removed': d_result['pct_removed'],
            }
            if cur_gid in s_gids_changed_domain:
                cur_result['classification'] = 'changed_domain'
            elif cur_gid in s_gids_no_markers:
                cur_result['classification'] = 'no_markers'
            elif cur_gid in s_gids_same_msa:
                cur_result['classification'] = 'same_msa'
            elif cur_gid in s_gids_run_on:
                cur_result['classification'] = 'run_on'
            else:
                raise Exception("???")
            if cur_gid in s_gids_run_on:
                if cur_gid in correct:
                    cur_result['tax_result'] = 'correct'
                elif cur_gid in wrong:
                    if cur_gid in gids_to_check:
                        cur_result['tax_result'] = 'check'
                    else:
                        cur_result['tax_result'] = 'congruent'
                else:
                    raise Exception('???')
            else:
                cur_result['tax_result'] = 'N/A'

            cur_result['markers_lost'] = ';'.join(d_result['markers_lost']),
            cur_result['markers_gained'] = ';'.join(d_result['markers_gained']),
            out.append(cur_result)

        df = pd.DataFrame(out)
        df.sort_values(by=['gid'], inplace=True, ignore_index=True)
        if not DEBUG:
            df.to_csv(self.output().path, sep='\t', index=False)
        return


def worker(job):
    gid, source, max_css, congruence, target_pct = job
    results = {
        'gid': gid,
        'source': source,
        'max_css': max_css
    }

    if source == 'gtdb':
        source = GuncRefDb.GTDB
    elif source == 'progenomes':
        source = GuncRefDb.PRO
    else:
        raise ValueError(f'Unknown source: {source}')

    genome = Genome(gid)
    contigs, pct_removed = genome.get_gunc_contigs_where_removed_equals_x_pct_genome_removed(
        pct=target_pct,
        max_congruence=congruence,
        max_css=max_css,
        source=source
    )
    results['pct_removed'] = pct_removed
    contigs = frozenset(contigs)

    gunc_tax = genome.get_gunc_max_css_inferred_taxonomy(max_css, source)
    results['gunc_tax'] = gunc_tax

    original_markers = genome.get_marker_hits()
    original_unq_hits = {**original_markers['muq'], **original_markers['unq']}
    original_unq_hits = set(original_unq_hits.keys())

    new_markers = genome.get_marker_hits(contigs)
    new_unq_hits = {**new_markers['muq'], **new_markers['unq']}
    new_unq_hits = set(new_unq_hits.keys())

    results['markers_lost'] = original_unq_hits - new_unq_hits
    results['markers_gained'] = new_unq_hits - original_unq_hits

    return results


def get_sister_taxa(tree: dendropy.Tree, gids: Set[str]) -> Dict[str, Set[str]]:
    out = dict()

    with tqdm(total=len({x.label for x in tree.taxon_namespace}.intersection(gids) )) as pbar:
        for leaf_node in tree.leaf_node_iter():
            if leaf_node.taxon.label in gids:
                out[leaf_node.taxon.label] = set()

                cur_node = leaf_node.parent_node
                while cur_node is not None:
                    support, taxon, auxiliary_info = parse_label(cur_node.label)

                    if taxon is not None and not taxon.startswith('s__'):
                        out[leaf_node.taxon.label].update({x.taxon.label for x in cur_node.leaf_nodes()})
                        pbar.update()
                        break
                    else:
                        cur_node = cur_node.parent_node

    return out



def match_nodes(untrimmed_tree: dendropy.Tree, trimmed_tree: dendropy.Tree, fail_gids: Set[str], d_gid_to_tax):

    # Calculate the descendant taxa to speed this process up
    log('Calculating descendant taxa (trimmed)')
    d_trimmed_desc_taxa = get_tree_node_to_desc_taxa(trimmed_tree)

    log('Calculating untrimmed descendant taxa')
    d_untrimmed_desc_taxa = get_tree_node_to_desc_taxa(untrimmed_tree)

    untrimmed_taxon_to_node = {x.taxon.label: x for x in untrimmed_tree.leaf_node_iter()}
    trimmed_taxon_to_node = {x.taxon.label: x for x in trimmed_tree.leaf_node_iter()}

    results = dict()

    for gid in tqdm(sorted(fail_gids), total=len(fail_gids)):

        # Get the node in the untrimmed tree
        untrimmed_node = untrimmed_taxon_to_node[gid]

        # Get the node in the trimmed tree
        trimmed_node = trimmed_taxon_to_node[f'TEST_{gid}']

        # Find the first parent node where no test genomes are present in any descendants
        cur_node = trimmed_node.parent_node
        trimmed_desc_taxa = set()
        trimmed_last_label_seen = None
        while cur_node is not None:
            support, taxon, auxiliary_info = parse_label(cur_node.label)

            if taxon is not None:
                trimmed_last_label_seen = taxon

            # desc_taxa = {x.taxon.label.replace('TEST_', '').replace('_C', '') for x in cur_node.leaf_nodes()}
            desc_taxa = {x for x in d_trimmed_desc_taxa[cur_node] if not x.startswith('TEST_')}
            n_desc_taxa = desc_taxa - {gid}
            if len(n_desc_taxa) > 0:
                trimmed_desc_taxa = desc_taxa
                break

            cur_node = cur_node.parent_node

        # How high up do we need to go in order to find the same descendant taxa in the untrimmed tree?
        n_steps = 0
        cur_node = untrimmed_node.parent_node
        untrimmed_last_label_seen = None
        while cur_node is not None:
            desc_taxa = d_untrimmed_desc_taxa[cur_node]

            support, taxon, auxiliary_info = parse_label(cur_node.label)

            if taxon is not None:
                untrimmed_last_label_seen = taxon

            if trimmed_desc_taxa.intersection(desc_taxa) == trimmed_desc_taxa:
                jaccard = len(trimmed_desc_taxa.intersection(desc_taxa)) / len(trimmed_desc_taxa.union(desc_taxa - fail_gids))
                tax_strings = [d_gid_to_tax[x] for x in desc_taxa]

                # find the highest rank that they agree at
                d_rank_counts = defaultdict(set)
                for tax_string in tax_strings:
                    ranks = tax_string.split(';')
                    for rank in ranks:
                        d_rank_counts[rank[0]].add(rank)

                highest_agree = 'd'
                for cur_rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
                    if len(d_rank_counts[cur_rank]) > 1:
                        break
                    highest_agree = cur_rank


                results[gid] = (n_steps, trimmed_last_label_seen, untrimmed_last_label_seen, jaccard, highest_agree, Counter(tax_strings))
                break

            cur_node = cur_node.parent_node
            n_steps += 1

    return results

def remove_species_from_tax_string(tax):
    out = list()
    for rank in tax.split(';'):
        if not rank.startswith('s__'):
            out.append(rank)
    return ';'.join(out)

def get_chimeric_halves_where_tax_strings_differ(d_trim, d_untrim, trim_tree, untrim_tree):

    d_trim_gid_to_node = {x.taxon.label: x for x in trim_tree.leaf_node_iter()}
    d_untrim_gid_to_node = {x.taxon.label: x for x in untrim_tree.leaf_node_iter()}

    fail_gids = {x.replace('TEST_', '').replace('_C', '') for x in d_trim if x.startswith('TEST_')}

    known_issues = {
        'GCA_900765805.1',
        'GCA_900759445.1',
        'GCA_905200565.1',
        'GCA_013372875.1',
        'GCA_012964215.1',
        'GCA_017525015.1',
        'GCA_014385135.2'
    }
    out = dict()

    gids_unsure = set()
    gids_same_sister = set()

    for gid in sorted(fail_gids):


        test_gid = f'TEST_{gid}'
        test_gid_c = f'TEST_{gid}_C'

        trim_tax_keep = d_trim[test_gid]
        # trim_tax_keep_no_s = remove_species_from_tax_string(trim_tax_keep)

        trim_tax_chim = d_trim[test_gid_c]
        # trim_tax_chim_no_s = remove_species_from_tax_string(trim_tax_chim)

        untrim_tax = d_untrim[gid]

        if trim_tax_keep != trim_tax_chim:

            test_gid_descs = get_descendants_of_first_named_node(d_trim_gid_to_node[test_gid])
            test_gid_descs_no_prefix = {x.replace('TEST_', '').replace('_C', '') for x in test_gid_descs}

            test_gid_sister = get_sister_taxon(d_trim_gid_to_node[test_gid])
            test_gid_sister_no_prefix = {x.replace('TEST_', '').replace('_C', '') for x in test_gid_sister}

            test_gid_c_descs = get_descendants_of_first_named_node(d_trim_gid_to_node[test_gid_c])
            test_gid_c_descs_no_prefix = {x.replace('TEST_', '').replace('_C', '') for x in test_gid_c_descs}

            test_gid_c_sister = get_sister_taxon(d_trim_gid_to_node[test_gid_c])
            test_gid_c_sister_no_prefix = {x.replace('TEST_', '').replace('_C', '') for x in test_gid_c_sister}

            untrim_descs = get_descendants_of_first_named_node(d_untrim_gid_to_node[gid])
            untrim_sister = get_sister_taxon(d_untrim_gid_to_node[gid])

            # If the sister taxa are the same between all three nodes it's unlikely they're a chimera
            if test_gid_sister_no_prefix == test_gid_c_sister_no_prefix == untrim_sister:
                gids_same_sister.add(gid)

            # If they're monophyletic within a specific taxon then it's unlikely they're chimeric, but may need investigation
            elif test_gid_descs_no_prefix == untrim_descs and test_gid_c_descs_no_prefix == untrim_descs:
                gids_unsure.add(gid)

            # Otherwise, manually inspect
            else:
                out[gid] = (trim_tax_keep, trim_tax_chim, untrim_tax)
    return out


def compare_tax_string_dicts(d_trim, d_untrim, trim_tree, untrim_tree):

    d_trim_gid_to_node = {x.taxon.label: x for x in trim_tree.leaf_node_iter()}
    d_untrim_gid_to_node = {x.taxon.label: x for x in untrim_tree.leaf_node_iter()}


    to_check = set()
    for test_gid, trim_tax in d_trim.items():

        if not test_gid.startswith('TEST_'):
            continue

        trim_tax_short, d_trim_tax, trim_genus = convert_phylorank_taxstring_to_dict(trim_tax)

        gid = test_gid.replace('TEST_', '').replace('_C', '')

        untrim_tax = d_untrim[gid]
        untrim_tax_short, d_untrim_tax, untrim_genus = convert_phylorank_taxstring_to_dict(untrim_tax)

        if trim_genus is not None and untrim_genus is not None:
            if trim_genus != untrim_genus:

                # Find all descendants of the first named node
                trim_descs = get_descendants_of_first_named_node(d_trim_gid_to_node[test_gid])
                untrim_descs = get_descendants_of_first_named_node(d_untrim_gid_to_node[gid])


                to_check.add(test_gid)
                print(test_gid)
                print(trim_tax)
                print(untrim_tax)
                print()
        else:
            to_check.add(test_gid)
            print('?')



    return


def convert_phylorank_taxstring_to_dict(tax):

    out = defaultdict(list)

    for rank in tax.split(';'):
        out[rank[0]].append(rank)

    lst = list()
    for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
        if len(out[rank]) == 1:
            lst.append(out[rank][0])
        else:
            break

    tax_genus = None
    if len(lst) >= 6:
        tax_genus = ';'.join(lst[:6])

    return ';'.join(lst), out, tax_genus


def get_sister_taxon(node: dendropy.Node):
    return {x.taxon.label for x in node.parent_node.leaf_nodes()}

def get_descendants_of_first_named_node(node: dendropy.Node):


    cur_node = node.parent_node
    while cur_node is not None:
        support, taxon, auxiliary_info = parse_label(cur_node.label)
        if taxon is not None:
            return {x.taxon.label for x in cur_node.leaf_nodes()}

        cur_node = cur_node.parent_node
    raise Exception(f'Unable to find node')



