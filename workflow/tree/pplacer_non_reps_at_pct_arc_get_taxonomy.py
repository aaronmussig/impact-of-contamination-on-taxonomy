import os

import dendropy
from luigi import LocalTarget
from tqdm import tqdm
import pandas as pd

from workflow.config import DEBUG, DIR_OUT_PPLACER_ARC_NON_REPS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.tree.pplacer_non_reps_at_pct_arc import PplacerNonRepsAtPctArc
from workflow.util.log import log

OUTPUT_DIRECTORY = DIR_OUT_PPLACER_ARC_NON_REPS


def get_taxonomy_from_tree(path_tree, df_meta):
    tree = dendropy.Tree.get_from_path(path_tree, schema='newick', preserve_underscores=True)

    taxa_to_check = {x.label for x in tree.taxon_namespace if x.label.startswith('GCA') or x.label.startswith('GCF')}
    print(f'{len(taxa_to_check):,} taxa to check')

    nodes_to_check = list()
    for leaf_node in tree.leaf_node_iter():
        if leaf_node.taxon.label in taxa_to_check:
            nodes_to_check.append(leaf_node)

    out = list()

    print(f'{len(nodes_to_check):,} leaf nodes to check')
    n_correct = 0
    n_wrong = 0

    for node in tqdm(nodes_to_check):
        sep_idx = node.taxon.label.rindex('_')
        gid, pct = node.taxon.label[:sep_idx], int(node.taxon.label[sep_idx + 1:])
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
            out.append({'gid': gid, 'pct': pct, 'tax': placed_tax})

            if placed_tax == expected_tax or placed_tax.split(';')[5] == expected_tax.split(';')[5]:
                n_correct += 1
            else:
                n_wrong += 1
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

            out.append({'gid': gid, 'pct': pct, 'tax': placed_tax})

            if placed_tax == expected_tax or congruent_rank == 5:
                n_correct += 1
            else:
                n_wrong += 1

    return out


class PplacerNonRepsAtPctArcGetTaxonomy(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'tree': PplacerNonRepsAtPctArc(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(OUTPUT_DIRECTORY, 'taxonomy.h5'))

    def run(self):
        log('Running pplacer for non arc representatives at pct values', title=True)
        self.make_output_dirs()

        self.make_output_dirs()
        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Getting taxonomy from tree')
        taxonomy = get_taxonomy_from_tree(self.input()['tree'].path, df_meta)

        log('Creating dataframe')
        df = pd.DataFrame(taxonomy)
        df.sort_values(by=['gid', 'pct'], inplace=True, ignore_index=True)

        if not DEBUG:
            self.save_hdf(df)
        return
