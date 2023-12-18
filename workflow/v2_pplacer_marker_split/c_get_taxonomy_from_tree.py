import os
import tempfile

import luigi
import pandas as pd
from magna.util.disk import move_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DIR_OUT_V2_PPLACER_MARKER_SPLIT, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetTsv
from workflow.util.log import log
from workflow.util.tree import read_taxonomy_from_tree, read_taxonomy_from_tree_pplacer
from workflow.v2_pplacer_marker_split.b_run_pplacer import V2PplacerMarkerSplitRunPplacer


class V2PplacerMarkerSplitGetTaxonomyFromTree(LuigiTask):
    remove_pct = luigi.FloatParameter(default=50)
    n_trees_output = luigi.IntParameter()

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_V2_PPLACER_MARKER_SPLIT, f'remove_pct_{self.remove_pct}')

    def requires(self):
        return {
            'pplacer': V2PplacerMarkerSplitRunPplacer(remove_pct=self.remove_pct, n_trees_output=self.n_trees_output),
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTargetTsv(os.path.join(self.root_dir, 'pplacer_taxonomy.tsv'))

    def run(self):
        log(f'V2PplacerMarkerSplitGetTaxonomyFromTree(remove_pct={self.remove_pct})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        non_rep_bac_gids = set(
            df_meta[(df_meta['gtdb_representative'] == 'f') & (df_meta['domain'] == 'd__Bacteria')].index)

        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        bac_non_rep_fail = non_rep_bac_gids.intersection(set(df_max_css.index))

        log('Processing trees')
        rows = list()
        for i in tqdm(range(self.n_trees_output), total=self.n_trees_output):
            tree = self.input()['pplacer'][str(i)].read()
            taxonomy = read_taxonomy_from_tree_pplacer(tree, d_gid_to_tax)

            for gid, tax in taxonomy.items():
                if gid.startswith('TEST_'):
                    rows.append({
                        'gid': gid[5:],
                        'tax': tax
                    })

        df = pd.DataFrame(rows)
        df.sort_values(by=['gid'], inplace=True, ignore_index=True)

        # Sanity check to make sure that all were actually captured
        if len(rows) != len(bac_non_rep_fail):
            raise Exception(f'Expected {len(bac_non_rep_fail)} non-representative bacteria, but only found {len(rows)}')

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp = os.path.join(tmp_dir, 'output.tsv')
            df.to_csv(path_tmp, sep='\t', index=False)

            log(f'Copying {get_file_size_fmt(path_tmp)} to {self.output().path}')
            if not DEBUG:
                move_file(path_tmp, self.output().path, checksum=True)

        return
