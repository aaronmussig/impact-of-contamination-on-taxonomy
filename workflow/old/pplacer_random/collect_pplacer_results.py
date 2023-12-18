import os
import sys

import dendropy
import pandas as pd
from tqdm import tqdm

from workflow.config import DIR_OUT_PPLACER_RANDOM_ARC, DIR_OUT_PPLACER_RANDOM_BAC, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.pplacer_random.generate_msas import GeneratePplacerRandomMsas
from workflow.util.log import log

sys.setrecursionlimit(10000)


def get_taxonomy_from_tree(tree, df_meta):
    taxa_to_check = {x.label for x in tree.taxon_namespace if len(x.label) > 18}
    log(f'{len(taxa_to_check):,} taxa to check')

    nodes_to_check = list()
    for leaf_node in tree.leaf_node_iter():
        if leaf_node.taxon.label in taxa_to_check:
            nodes_to_check.append(leaf_node)

    out = list()

    log(f'{len(nodes_to_check):,} leaf nodes to check')
    n_correct = 0
    n_wrong = 0
    lst_wrong = list()

    for node in tqdm(nodes_to_check):
        rep_id, gid_id, gid, pct = node.taxon.label.split('__')
        rep_id, gid_id, pct = int(rep_id), int(gid_id), int(pct)
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

        if len(closest_gids) == 1:
            placed_rep = closest_gids[0]
            placed_tax = df_meta.loc[placed_rep]['gtdb_taxonomy']
            out.append({'gid': gid, 'rep': rep_id, 'gid_id': gid_id, 'pct': pct, 'tax': placed_tax})

            if placed_tax == expected_tax or placed_tax.split(';')[5] == expected_tax.split(';')[5]:
                n_correct += 1
            else:
                n_wrong += 1
                lst_wrong.append(f'{gid} {expected_tax} {placed_tax}')
        else:

            # find the most congruent taxonomy
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

            out.append({'gid': gid, 'rep': rep_id, 'gid_id': gid_id, 'pct': pct, 'tax': placed_tax})

            if placed_tax == expected_tax or congruent_rank == 5:
                n_correct += 1
            else:
                n_wrong += 1
                lst_wrong.append(f'{gid} {expected_tax} {placed_tax}')

    log(f'{n_correct:,} correct, {n_wrong:,} wrong')
    if n_wrong > 0:
        [log(x) for x in lst_wrong]
    return out


def validate_all_gids_placed(path_msa, tree_taxa):
    msa_gids = set()
    with open(path_msa) as f:
        for line in f.readlines():
            if not line.startswith('>'):
                continue
            msa_gids.add(line.strip()[1:])

    if msa_gids != tree_taxa:
        raise Exception('Missing gids!')
    return


def get_results_from_root_directory(dir_name, df_meta):
    out = list()
    for batch in os.listdir(dir_name):
        path_msa = os.path.join(dir_name, batch, 'input.fasta')
        path_tree = os.path.join(dir_name, batch, 'pplacer_random.tree')

        if not os.path.isfile(path_tree) and DEBUG:
            continue

        tree = dendropy.Tree.get_from_path(path_tree, schema='newick', preserve_underscores=True)
        placed_taxa = {x.label for x in tree.taxon_namespace if len(x.label) > 18}

        validate_all_gids_placed(path_msa, placed_taxa)

        tax_from_tree = get_taxonomy_from_tree(tree, df_meta)

        out.append(tax_from_tree)
    log(f'Enqueued {len(out):,} msas from {dir_name}')
    return out


class CollectPplacerRandomResults(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            '_pplacer_msas': GeneratePplacerRandomMsas(),
        }

    def output(self):
        return {
            'arc': LocalTargetHdf5(os.path.join(DIR_OUT_PPLACER_RANDOM_ARC, 'arc_results.h5')),
            'bac': LocalTargetHdf5(os.path.join(DIR_OUT_PPLACER_RANDOM_BAC, 'bac_results.h5')),
        }

    def run(self):
        log('Collecting pplacer results (random removal)', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        arc_input = get_results_from_root_directory(DIR_OUT_PPLACER_RANDOM_ARC, df_meta)
        arc_rows = list()
        [arc_rows.extend(x) for x in arc_input]
        arc_df = pd.DataFrame(arc_rows)
        arc_df.sort_values(by=['rep', 'gid_id'], inplace=True)

        bac_input = get_results_from_root_directory(DIR_OUT_PPLACER_RANDOM_BAC, df_meta)
        bac_rows = list()
        [bac_rows.extend(x) for x in bac_input]
        bac_df = pd.DataFrame(bac_rows)
        bac_df.sort_values(by=['rep', 'gid_id'], inplace=True)

        log('Saving dataframes')
        if not DEBUG:
            self.save_hdf(arc_df, path=self.output()['arc'].path)
            self.save_hdf(bac_df, path=self.output()['bac'].path)

        return
