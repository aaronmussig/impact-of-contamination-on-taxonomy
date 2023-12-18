import os

import dendropy
import luigi
import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, FASTTREE_MARKER_ANALYSIS_N_BATCHES
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_tree import GtdbTreeBacR207
from workflow.fasttree_marker_analysis.random_d_analyse_decorated import FastTreeMarkerAnalyseDecoratedRandom
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.tree import parse_label


class FastTreeMarkerAnalyseDecoratedRandomGetParentBs(LuigiTask):
    target_pct = luigi.FloatParameter()

    def requires(self):
        out = {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'r207_bac_ref_tree': GtdbTreeBacR207(),
        }

        for i in range(FASTTREE_MARKER_ANALYSIS_N_BATCHES):
            out[f'analysis_{i}'] = FastTreeMarkerAnalyseDecoratedRandom(batch_id=i, target_pct=self.target_pct)
        return out

    def output(self):
        return LocalTargetHdf5(
            os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, 'random_e_analyse_decorated_get_parent_bs.h5'))

    def run(self):
        log(f'Analysing results of decorated tree (random) (named node)  (pct={self.target_pct})',
            title=True)
        # self.make_output_dirs()

        log('Loading BAC120 R207 reference tree')
        r207_bac_ref_tree = self.input()['r207_bac_ref_tree'].read()
        for leaf_node in r207_bac_ref_tree.leaf_node_iter():
            leaf_node.taxon.label = leaf_node.taxon.label[3:]

        log('Loading previous results')
        incongruent_values = list()
        congruent_values = list()
        no_markers_values = list()

        for batch_id in range(3 if DEBUG else FASTTREE_MARKER_ANALYSIS_N_BATCHES):
            log(f'Processing batch id={batch_id}')
            df = pd.read_csv(self.input()[f'analysis_{batch_id}'].path, sep='\t')
            df.set_index('gid', inplace=True)

            # log('Breaking into sets')
            gids_congruent = set()
            gids_incongruent = set()
            gids_no_markers = set()

            # Congruent
            gids_congruent.update(set(df[df['classification'] == 'same_msa'].index))
            gids_congruent.update(set(df[(df['classification'] == 'run_on') & (df['tax_result'] == 'correct')].index))
            gids_congruent.update(set(df[(df['classification'] == 'run_on') & (df['tax_result'] == 'congruent')].index))

            # Incongruent
            gids_incongruent.update(set(df[(df['classification'] == 'run_on') & (df['tax_result'] == 'check')].index))
            gids_incongruent.update(set(df[df['classification'] == 'changed_domain'].index))

            # No markers
            gids_no_markers.update(set(df[df['classification'] == 'no_markers'].index))

            # log('Sanity check')
            assert (gids_congruent.union(gids_no_markers).union(gids_incongruent) == set(df.index))
            assert (len(gids_congruent.intersection(gids_no_markers).intersection(gids_incongruent)) == 0)

            # log('Checking bootstrap values in reference tree')
            incongruent_values.extend(check_ref_tree(r207_bac_ref_tree, gids_incongruent, 'incongruent', batch_id))
            congruent_values.extend(check_ref_tree(r207_bac_ref_tree, gids_congruent, 'congruent', batch_id))
            no_markers_values.extend(check_ref_tree(r207_bac_ref_tree, gids_no_markers, 'no_markers', batch_id))

        log('Creating dataframe')
        df = pd.DataFrame(incongruent_values + congruent_values + no_markers_values)
        df.sort_values(by=['batch_id', 'set', 'gid', 'bs'], inplace=True, ignore_index=True)

        if not DEBUG:
            self.save_hdf(df)

        return


def check_ref_tree(ref_tree: dendropy.Tree, gid_set, set_name, batch_id):
    out = list()

    for gid in tqdm(sorted(gid_set)):
        highest_taxon = None
        highest_bs = None

        leaf = ref_tree.find_node_with_taxon_label(gid)
        cur_parent = leaf.parent_node
        while cur_parent is not None:

            support, taxon, auxiliary_info = parse_label(cur_parent.label)

            if taxon is not None:
                highest_bs = support
                highest_taxon = taxon
                break

            cur_parent = cur_parent.parent_node

        if highest_bs is None or highest_taxon is None:
            out.append({
                'gid': gid,
                'batch_id': batch_id,
                'bs': -100,
                'set': set_name,
                'taxon': 'N/A',
            })

        else:
            out.append({
                'batch_id': batch_id,
                'set': set_name,
                'gid': gid,
                'bs': int(highest_bs),
                'taxon': highest_taxon,
            })

        if DEBUG and len(out) > 20:
            break

    return out
