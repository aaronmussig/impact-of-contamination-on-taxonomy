import multiprocessing as mp
import os
from typing import Set, Dict

import dendropy
import luigi
import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_tree import GtdbTreeBacR207
from workflow.fasttree_full_tree_non_reps.gunc_c_decorate_fasttree import \
    FastTreeFullSetNonRepsCreateBatchesGuncDecorate
from workflow.fasttree_marker_analysis.random_d_analyse_decorated import decorated_tax_to_agreed, get_taxon_to_gids
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


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


class FastTreeFullSetNonRepsCreateBatchesGuncDecorateAnalysis(LuigiTask):
    congruence = luigi.FloatParameter()
    target_pct = luigi.FloatParameter()

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, f'c{str(self.congruence)}_p{str(self.target_pct)}')

    def requires(self):
        return {
            'decorated': FastTreeFullSetNonRepsCreateBatchesGuncDecorate(congruence=self.congruence,
                                                                         target_pct=self.target_pct),
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'r207_bac_ref_tree': GtdbTreeBacR207(),
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, 'decorated_analysis_results2.tsv'))

    def run(self):
        log(f'FastTreeFullSetNonRepsCreateBatchesGuncDecorateAnalysis(congruence={self.congruence}) (pct={self.target_pct})',
            title=True)
        self.make_output_dirs()

        log('Loading BAC120 R207 reference tree')
        r207_bac_ref_tree = self.input()['r207_bac_ref_tree'].read()
        for leaf_node in r207_bac_ref_tree.leaf_node_iter():
            leaf_node.taxon.label = leaf_node.taxon.label[3:]

        log('Loading decorated tree')
        decorated_tree = self.input()['decorated'].read()

        # Set paths
        phylo_dir = os.path.dirname(self.input()['decorated'].path)
        root_dir = os.path.dirname(phylo_dir)

        log('Loading gids that changed domain')
        s_gids_changed_domain = read_gids_file(os.path.join(root_dir, 'gids_changed_domain.txt'))
        s_gids_no_markers = read_gids_file(os.path.join(root_dir, 'gids_no_markers.txt'))
        s_gids_same_msa = read_gids_file(os.path.join(root_dir, 'gids_same_msa.txt'))
        s_gids_run_on = read_gids_file(os.path.join(root_dir, 'gids_to_run_on.txt'))
        s_all_gids = s_gids_changed_domain.union(s_gids_no_markers).union(s_gids_same_msa).union(s_gids_run_on)
        total_gids = len(s_all_gids)
        log(f'Loaded {total_gids:,} gids')

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_taxonomy = df_meta['gtdb_taxonomy'].to_dict()

        path_decorated_tax = os.path.join(phylo_dir, 'decorated.tree-taxonomy')
        d_gid_to_decorated_tax = read_tax_file(path_decorated_tax)

        log('Getting taxonomic novelty (reps only)')
        d_taxon_to_gids = get_taxon_to_gids(df_meta[df_meta['gtdb_representative'] == 't']['gtdb_taxonomy'].to_dict())

        correct = set()
        wrong = set()
        for test_gid, decorated_tax in tqdm(sorted(d_gid_to_decorated_tax.items())):
            if test_gid.startswith('TEST_'):
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
                    wrong.add(gid)
                else:
                    if not true_tax.endswith(d_decorated_tax_suffix):
                        wrong.add(gid)
                    else:
                        correct.add(gid)

                # print('mismatch', test_gid, true_tax, decorated_tax)
            else:
                correct.add(gid)

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
    for leaf_node in tree.leaf_node_iter():
        if leaf_node.taxon.label in gids:
            out[leaf_node.taxon.label] = set()
            for sister_node in leaf_node.parent_node.leaf_nodes():
                out[leaf_node.taxon.label].add(sister_node.taxon.label)
    assert (len(out) == len(gids))
    return out
