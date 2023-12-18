import os
from collections import defaultdict

import luigi
import pandas as pd

from workflow.config import DEBUG, DIR_OUT_V2_FASTTREE_MARKER_SPLIT
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask, LocalTargetTsv
from workflow.util.log import log
from workflow.v2_fasttree_marker_split.d_decorate_fasttree import V2FastTreeMarkerSplitDecorateFastTree


class V2FastTreeMarkerSplitAnalyseDecoratedFastTree(LuigiTask):
    remove_pct = luigi.FloatParameter(default=50)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_V2_FASTTREE_MARKER_SPLIT, f'remove_pct_{self.remove_pct}')

    def requires(self):
        return {
            'decorated': V2FastTreeMarkerSplitDecorateFastTree(remove_pct=self.remove_pct),
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTargetTsv(os.path.join(self.root_dir, 'decorated_analysis_results.tsv'))

    def run(self):
        log(f'V2FastTreeMarkerSpliAnalyseDecoratedFastTree(remove_pct={self.remove_pct})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()

        log('Loading decorated trimmed tree')
        tree_trimmed = self.input()['decorated'].read()

        log('Loading tax strings from decorated trimmed tree')
        # Note: If there are duplicate ranks, then the first one that appears was closest to the leaf node
        # d_test_gid_to_tax = read_taxonomy_from_tree(tree_trimmed)
        trim_root_dir = os.path.dirname(self.input()['decorated'].path)
        trim_decorated_tax = os.path.join(trim_root_dir, 'decorated.tree-taxonomy')

        d_test_gid_to_tax = read_taxonomy_from_tax_file(trim_decorated_tax)
        df = compare_tax_strings_from_tree(d_test_gid_to_tax, d_gid_to_tax)

        df.sort_values(by=['gid'], inplace=True, ignore_index=True)
        if not DEBUG:
            df.to_csv(self.output().path, sep='\t', index=False)
        return


def read_taxonomy_from_tax_file(path):
    d_gid_to_decorated_tax = dict()
    with open(path) as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            cur_tax = tax.replace('; ', ';')
            d_cur_tax = parse_phylorank_tax_string(cur_tax)
            for rank in 'dpcofgs':
                if rank not in d_cur_tax:
                    d_cur_tax[rank] = [f'{rank}__']
            cur_out = list()
            for rank in 'dpcofgs':
                cur_out.extend(d_cur_tax[rank])
            d_gid_to_decorated_tax[gid] = ';'.join(cur_out)
    return d_gid_to_decorated_tax


def compare_tax_strings_from_tree(d_gid_to_tax_test, d_gid_to_tax_true):
    out = list()

    for gid, tax_test in d_gid_to_tax_test.items():
        if not gid.startswith('TEST_'):
            continue

        tax_true = d_gid_to_tax_true[gid[5:]]
        d_tax_true = {x[0]: x for x in tax_true.split(';')}

        d_cur_out = {
            'gid': gid,
            'test_tax': tax_test,
            'true_tax': tax_true,
            'highest_disagree_rank': 'na',
            'disagree_type': 'na',
            'agree': True,
        }

        if tax_test != tax_true:

            # In-case phylorank has nested/skipped ranks
            test_prefixes = [x[0] for x in tax_test.split(';')]
            if not (len(test_prefixes) == 7 and len(set(test_prefixes)) == 7):
                d_test_tax = parse_phylorank_tax_string(tax_test)

                for rank in 'dpcofgs':
                    # literal compare
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

                    if test_rank != true_rank:
                        d_cur_out['agree'] = False
                        d_cur_out['highest_disagree_rank'] = test_rank[0]

                        if len(test_rank) == 3:
                            d_cur_out['disagree_type'] = 'missing'
                        else:
                            d_cur_out['disagree_type'] = 'mismatch'
                        break


            # Simple case where you can just compare the tax strings literally to find mismatches
            else:
                for test_rank, true_rank in zip(tax_test.split(';'), tax_true.split(';')):
                    if test_rank.startswith('s__'):
                        continue

                    if test_rank != true_rank:
                        d_cur_out['agree'] = False
                        d_cur_out['highest_disagree_rank'] = test_rank[0]

                        if len(test_rank) == 3:
                            d_cur_out['disagree_type'] = 'missing'
                        else:
                            d_cur_out['disagree_type'] = 'mismatch'
                        break

        # Save the row
        out.append(d_cur_out)

    df = pd.DataFrame(out)

    return df


def parse_phylorank_tax_string(taxonomy):
    out = defaultdict(list)
    for rank in taxonomy.split(';'):
        out[rank[0]].append(rank)
    for rank in 'dpcofgs':
        if rank not in out:
            out[rank] = [f'{rank}__']
    return out
