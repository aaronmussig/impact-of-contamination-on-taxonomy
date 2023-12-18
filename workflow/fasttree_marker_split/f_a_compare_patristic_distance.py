import os
from collections import defaultdict

import dendropy
import luigi
import numpy as np
import pandas as pd
from phylodm import PhyloDM
from scipy import stats
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTTREE_MARKER_SPLIT, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_full_tree_and_failed_no_trim.b_run_fasttree import FastTreeFullTreeAndFailedNoTrimRunFastTree
from workflow.fasttree_marker_split.b_run_fasttree import FastTreeMarkerSplitRunFastTree
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class FastTreeMarkerSplitComparePatristicDistance(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_MARKER_SPLIT, f'pct_{self.target_pct}')

    def requires(self):
        return {
            'split_tree': FastTreeMarkerSplitRunFastTree(target_pct=self.target_pct),
            'ref_tree': FastTreeFullTreeAndFailedNoTrimRunFastTree(replicate_id=0),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, 'patristic_distances.h5'))

    def run(self):
        log(f'FastTreeMarkerSplitComparePatristicDistance(p={self.target_pct})', title=True)
        self.make_output_dirs()

        log('Calculating PDM for split tree')
        split_pdm = generate_pdm(self.input()['split_tree'].path)
        split_taxon_to_idx = {t: i for i, t in enumerate(split_pdm.taxa())}
        fail_gids = frozenset(
            {x.replace('TEST_', '') for x in split_taxon_to_idx if x.startswith('TEST_') and not x.endswith('_C')})

        log('Loading reference PDM')
        ref_pdm = generate_pdm(self.input()['ref_tree'].path)

        nn = 10000
        log(f'Getting nearest {nn} nearest neighbours for each fail genome in the reference tree')
        nearest_neighbours = get_nearest_neighbours(fail_gids, nn, ref_pdm)

        log('Calculating nearest neighbour distances')
        d_distances = compare_nearest_neighbour_distances(nearest_neighbours, fail_gids, ref_pdm, split_pdm)

        log('Calculating KS statistics')
        d_ks_results = ks_test(d_distances)

        log('Calculating distance to each genomes chimeric half')
        split_chim_deltas = calculate_delta_to_chimeric_half(split_pdm, fail_gids)

        log('Merging results')
        rows = list()
        for gid in fail_gids:
            d_ks_result = d_ks_results[gid]
            split_chim_delta = split_chim_deltas[gid]
            rows.append({
                'gid': gid,
                'ks_ref_vs_keep': d_ks_result['ref_vs_keep'],
                'ks_ref_vs_chim': d_ks_result['ref_vs_chim'],
                'ks_keep_vs_chim': d_ks_result['keep_vs_chim'],
                'pd_between_halves': split_chim_delta
            })
        df = pd.DataFrame(rows)

        if not DEBUG:
            log('Saving dataframe')
            self.save_hdf(df)

        return


def generate_pdm(path_tree: str):
    log('Loading tree')
    tree = dendropy.Tree.get_from_path(path_tree, schema='newick', preserve_underscores=True)
    n_updated = 0
    for edge in tree.postorder_edge_iter():
        if edge.length is None and not tree.seed_node.edge == edge:  # TODO: Make sure this is not the seed node edge
            edge.length = 0.0
            n_updated += 1
    log('Updated {} edges with missing lengths'.format(n_updated))
    log('Creating PDM')
    pdm = PhyloDM.load_from_dendropy(tree)

    log('Creating row vector')
    pdm.compute_row_vec()
    return pdm


def calculate_delta_to_chimeric_half(pdm, fail_gids):
    out = dict()
    norm = pdm.length()
    for gid in tqdm(sorted(fail_gids)):
        out[gid] = pdm.distance(f'TEST_{gid}', f'TEST_{gid}_C') / norm
    return out


def get_nearest_neighbours(taxa, n, pdm):
    out = dict()

    dm = pdm.dm(norm=True)
    idx_to_taxon = pdm.taxa()
    taxon_to_idx = {t: i for i, t in enumerate(idx_to_taxon)}

    for taxon in sorted(taxa):
        cur_idx = taxon_to_idx[taxon]
        closet_idxs = np.argsort(dm[cur_idx, :])[:n + 1]
        closest_gids = [idx_to_taxon[i] for i in closet_idxs if i != cur_idx][:n]
        out[taxon] = closest_gids

    return out


def compare_nearest_neighbour_distances(nearest_neighbours, fail_gids, ref_pdm, split_pdm):
    split_norm = split_pdm.length()
    ref_norm = ref_pdm.length()

    out = defaultdict(list)

    for taxon, lst_nn in tqdm(nearest_neighbours.items()):
        n_added = 0

        for nn_idx, nn in enumerate(lst_nn):
            if n_added >= 100:
                break

            if nn not in fail_gids:
                # Get the distance in the reference tree
                ref_dist = ref_pdm.distance(taxon, nn) / ref_norm

                # Get the distance in the query tree to the keep half
                keep_dist = split_pdm.distance(f'TEST_{taxon}', nn) / split_norm
                chim_dist = split_pdm.distance(f'TEST_{taxon}_C', nn) / split_norm

                out[taxon].append({
                    'taxon': taxon,
                    'nn': nn,
                    'nn_idx': nn_idx,
                    'ref_dist': ref_dist,
                    'keep_dist': keep_dist,
                    'chim_dist': chim_dist,
                })
                n_added += 1

    return out


def ks_test(d_distances):
    out = dict()
    for taxon, rows in tqdm(d_distances.items()):
        df_rows = pd.DataFrame(rows)

        ks_ref_v_keep = stats.kstest(df_rows['ref_dist'].values, df_rows['keep_dist'].values)
        ks_ref_v_chim = stats.kstest(df_rows['ref_dist'].values, df_rows['chim_dist'].values)
        ks_keep_v_chim = stats.kstest(df_rows['keep_dist'].values, df_rows['chim_dist'].values)

        out[taxon] = {
            'ref_vs_keep': ks_ref_v_keep.pvalue,
            'ref_vs_chim': ks_ref_v_chim.pvalue,
            'keep_vs_chim': ks_keep_v_chim.pvalue,
        }

    return out
