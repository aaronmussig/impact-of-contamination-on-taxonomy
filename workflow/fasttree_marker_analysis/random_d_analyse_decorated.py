import multiprocessing as mp
import os
from typing import Set, Dict, List

import dendropy
import luigi
import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_tree import GtdbTreeBacR207
from workflow.fasttree_marker_analysis.random_c_decorate_fasttree import FastTreeMarkerAnalysisDecorateRandomFastTree
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
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


class FastTreeMarkerAnalyseDecoratedRandom(LuigiTask):
    target_pct = luigi.FloatParameter()
    batch_id = luigi.IntParameter()

    def requires(self):
        return {
            'decorated': FastTreeMarkerAnalysisDecorateRandomFastTree(batch_id=self.batch_id,
                                                                      target_pct=self.target_pct),
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'r207_bac_ref_tree': GtdbTreeBacR207(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, f'{self.batch_id}_{self.target_pct}',
                                        'decorated_analysis_results_2.tsv'))

    def run(self):
        log(f'Analysing results of decorated tree (random) (batch_id={self.batch_id}) (pct={self.target_pct})',
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
        path_contigs_removed = os.path.join(root_dir, 'gids_contigs_removed.txt')

        d_contigs_removed = read_contigs_removed(path_contigs_removed)

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

        log('Getting taxonomic novelty (reps only)')
        d_taxon_to_gids = get_taxon_to_gids(df_meta[df_meta['gtdb_representative'] == 't']['gtdb_taxonomy'].to_dict())

        path_decorated_tax = os.path.join(phylo_dir, 'decorated.tree-taxonomy')
        d_gid_to_decorated_tax = read_tax_file(path_decorated_tax)
        d_gid_to_decorated_tax_strip = {k.replace('TEST_', ''): v for k, v in d_gid_to_decorated_tax.items()}

        correct = set()
        wrong = set()
        for test_gid, decorated_tax in sorted(d_gid_to_decorated_tax.items()):
            if test_gid.startswith('TEST_'):
                gid = test_gid[5:]
            else:
                continue
            true_tax = d_gid_to_taxonomy[gid]
            if true_tax != decorated_tax:

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

        # Get the equivalent node in the true reference tree
        ref_sister_taxa = get_sister_taxa(r207_bac_ref_tree, wrong)
        dec_sister_taxa = get_sister_taxa(decorated_tree, {f'TEST_{x}' for x in wrong})

        # Compare the results
        gids_to_check = set()
        for gid in wrong:
            ref_taxa = ref_sister_taxa[gid]
            dec_taxa = dec_sister_taxa[f'TEST_{gid}']
            dec_taxa_strip = {x.replace('TEST_', '') for x in dec_taxa}

            if dec_taxa_strip != ref_taxa:
                gids_to_check.add(gid)

        log('Creating queue')
        queue = list()
        for gid, (pct_removed, contigs_removed) in d_contigs_removed.items():
            queue.append((gid, contigs_removed, pct_removed))

        log('Processing queue')
        with mp.Pool(processes=int(mp.cpu_count() * 0.8)) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        out = list()
        for d_result in results:
            cur_gid = d_result['gid']
            cur_result = {
                'gid': cur_gid,
                'batch_id': self.batch_id,
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

        log(f'TO CHECK: {len(gids_to_check)}')
        log(f'Saving to: {self.output().path}')

        df = pd.DataFrame(out)
        df.sort_values(by=['gid'], inplace=True, ignore_index=True)
        if not DEBUG:
            df.to_csv(self.output().path, sep='\t', index=False)
        return


def worker(job):
    gid, s_contigs_removed, pct_removed = job

    results = {
        'gid': gid,
    }

    genome = Genome(gid)

    results['pct_removed'] = pct_removed

    if pct_removed > 0:

        original_markers = genome.get_marker_hits()
        original_unq_hits = {**original_markers['muq'], **original_markers['unq']}
        original_unq_hits = set(original_unq_hits.keys())

        new_markers = genome.get_marker_hits(s_contigs_removed)
        new_unq_hits = {**new_markers['muq'], **new_markers['unq']}
        new_unq_hits = set(new_unq_hits.keys())

        results['markers_lost'] = original_unq_hits - new_unq_hits
        results['markers_gained'] = new_unq_hits - original_unq_hits

    else:
        results['markers_lost'] = set()
        results['markers_gained'] = set()

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


def read_contigs_removed(path):
    out = dict()
    with open(path, 'r') as f:
        for line in f.readlines():
            cols = line.strip().split('\t')
            if len(cols) == 3:
                gid, pct_removed, contigs_removed = cols
                s_contigs_removed = contigs_removed.split(';')
                s_contigs_removed = frozenset(s_contigs_removed)
            else:
                gid, pct_removed = cols
                s_contigs_removed = frozenset()
            out[gid] = (float(pct_removed), s_contigs_removed)
    return out


def get_taxon_to_gids(d_gid_to_taxonomy) -> Dict[str, Set[str]]:
    out = dict()
    for gid, taxonomy in d_gid_to_taxonomy.items():
        for taxon in taxonomy.split(';'):
            if taxon not in out:
                out[taxon] = set()
            out[taxon].add(gid)
    return out


def decorated_tax_to_dict(tax: str) -> Dict[str, List[str]]:
    out = dict()
    for taxon in tax.split(';'):
        rank = taxon[0]
        if rank not in out:
            out[rank] = list()
        out[rank].append(taxon)
    return out


def decorated_tax_to_agreed(tax: str):
    d_tax = decorated_tax_to_dict(tax)
    out = list()
    for rank in reversed('dpcofgs'):
        tax_cnt = len(d_tax.get(rank, list()))
        if tax_cnt != 1:
            break
        out.append(d_tax[rank][0])
    return ';'.join(reversed(out))
