import os
import sys
import tempfile

import luigi
import pandas as pd
from luigi import LocalTarget
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DIR_OUT_PPLACER_MARKER_ANALYSIS, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask, LocalTargetTsv
from workflow.pplacer_marker_analysis.a_create_batches import PplacerMarkerAnalysisCreateBatches
from workflow.pplacer_marker_analysis.b_run_pplacer import PplacerMarkerAnalysisRunPplacer
from workflow.util.log import log

sys.setrecursionlimit(10000)


class PplacerMarkerAnalysisAnalysePplacerResults(LuigiTask):
    # Required for both
    domain = luigi.ChoiceParameter(choices=['arc', 'bac'])
    target_pct = luigi.FloatParameter()

    # Required for GUNC-informed
    congruence = luigi.FloatParameter(default=-1)

    # Required for random
    random = luigi.BoolParameter(default=False)
    batch_id = luigi.IntParameter(default=-1)

    def requires(self):
        return {
            'msa': PplacerMarkerAnalysisCreateBatches(
                domain=self.domain,
                target_pct=self.target_pct,
                congruence=self.congruence,
                random=self.random,
                batch_id=self.batch_id
            ),
            'tree': PplacerMarkerAnalysisRunPplacer(
                domain=self.domain,
                target_pct=self.target_pct,
                congruence=self.congruence,
                random=self.random,
                batch_id=self.batch_id
            ),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        if self.random:
            target_dir = os.path.join(
                DIR_OUT_PPLACER_MARKER_ANALYSIS,
                f'random_{self.domain}_p{self.target_pct}_b{self.batch_id}'
            )
        else:
            target_dir = os.path.join(
                DIR_OUT_PPLACER_MARKER_ANALYSIS,
                f'gunc_{self.domain}_p{self.target_pct}_c{self.congruence}'
            )
        return LocalTargetTsv(os.path.join(target_dir, 'results.tsv'))

    def run(self):
        if self.random:
            log(f'RANDOM: PplacerMarkerAnalysisAnalysePplacerResults('
                f'pct={self.target_pct}, batch_id={self.batch_id}, '
                f'domain={self.domain})', title=True)
        else:
            log(f'GUNC: PplacerMarkerAnalysisAnalysePplacerResults('
                f'congruence={self.congruence}, pct={self.target_pct}, '
                f'domain={self.domain})', title=True)
        self.make_output_dirs()

        df_msa_log = self.input()['msa']['log'].maybe_read_cached()
        pplacer_tree = self.input()['tree'].read()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Getting taxonomy from tree')
        rows = get_taxonomy_from_tree(pplacer_tree, df_meta)

        df_tax = pd.DataFrame(rows)
        df_tax.set_index('gid', inplace=True)

        df_merged = df_msa_log.merge(df_tax, how='left', left_index=True, right_index=True)
        df_merged.reset_index(inplace=True)
        df_merged.rename(columns={"index": "gid"}, inplace=True)
        df_merged.sort_values(by='gid', inplace=True, ignore_index=True)

        if len(df_merged) != len(df_msa_log):
            raise Exception('?')

        if not DEBUG:
            log(f'Saving to: {self.output().path}')
            with tempfile.TemporaryDirectory() as tmp_dir:
                path_tmp = os.path.join(tmp_dir, 'results.tsv')
                df_merged.to_csv(path_tmp, sep='\t', index=False)
                copy_file(path_tmp, self.output().path, checksum=True)

        return


def get_taxonomy_from_tree(tree, df_meta):
    taxa_to_check = {x.label for x in tree.taxon_namespace if x.label.startswith("TEST_")}
    log(f'{len(taxa_to_check):,} taxa to check')

    nodes_to_check = dict()
    for leaf_node in tree.leaf_node_iter():
        if leaf_node.taxon.label in taxa_to_check:
            nodes_to_check[leaf_node.taxon.label[5:]] = leaf_node

    out = list()

    log(f'{len(nodes_to_check):,} leaf nodes to check')
    n_correct = 0
    n_wrong = 0
    lst_wrong = list()

    for gid, node in tqdm(nodes_to_check.items()):
        expected_tax = df_meta.loc[gid]['gtdb_taxonomy']

        # take the full taxonomy of this genome to be the same
        parent_node = node.parent_node
        closest_gids = list()
        while parent_node is not None:
            closest_gids = [x.taxon.label[3:] for x in parent_node.leaf_nodes() if
                            x.taxon.label.startswith('GB_') or x.taxon.label.startswith('RS_')]
            if len(closest_gids) == 0:
                parent_node = parent_node.parent_node
            else:
                break

        # Easy case - just take the sister taxonomy
        if len(closest_gids) == 1:
            placed_rep = closest_gids[0]
            placed_tax = df_meta.loc[placed_rep]['gtdb_taxonomy']

            if placed_tax == expected_tax or placed_tax.split(';')[5] == expected_tax.split(';')[5]:
                n_correct += 1
                is_correct = True
            else:
                n_wrong += 1
                is_correct = False
                lst_wrong.append(f'{gid} {expected_tax} {placed_tax}')

            out.append({
                'gid': gid,
                'tax': placed_tax,
                'correct': is_correct
            })

        # find the most congruent taxonomy
        else:
            taxonomies = [df_meta.loc[x]['gtdb_taxonomy'].split(';') for x in closest_gids]
            congruent_rank = 0
            for i in range(0, 7):
                taxa = {x[i] for x in taxonomies}
                if len(taxa) > 1:
                    break
                congruent_rank = i

            new_tax = list()
            for i, rank in enumerate(taxonomies[0]):
                if i > congruent_rank:
                    new_tax.append(rank[0:3])
                else:
                    new_tax.append(rank)
            placed_tax = ';'.join(new_tax)

            if placed_tax == expected_tax or congruent_rank == 5:
                n_correct += 1
                is_correct = True
            else:
                n_wrong += 1
                is_correct = False
                lst_wrong.append(f'{gid} {expected_tax} {placed_tax}')

            out.append({
                'gid': gid,
                'tax': placed_tax,
                'correct': is_correct
            })

    log(f'{n_correct:,} correct, {n_wrong:,} wrong')
    if n_wrong > 0:
        [log(x) for x in lst_wrong]
    return out
