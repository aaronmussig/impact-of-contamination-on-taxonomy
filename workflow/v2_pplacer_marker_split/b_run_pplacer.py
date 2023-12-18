import os
import tempfile

import luigi
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DIR_OUT_V2_PPLACER_MARKER_SPLIT
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.model.pplacer import Pplacer
from workflow.util.log import log
from workflow.v2_pplacer_marker_split.a_create_msa import V2PplacerMarkerSplitCreateMsa


class V2PplacerMarkerSplitRunPplacer(LuigiTask):
    remove_pct = luigi.FloatParameter(default=50)
    n_trees_output = luigi.IntParameter()

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_V2_PPLACER_MARKER_SPLIT, f'remove_pct_{self.remove_pct}')

    def requires(self):
        return {
            'msa': V2PplacerMarkerSplitCreateMsa(remove_pct=self.remove_pct)
        }

    def output(self):
        return {str(x): LocalTargetTree(os.path.join(self.root_dir, f'pplacer_{x}.tree')) for x in
                range(self.n_trees_output)}

    def run(self):
        log(f'V2PplacerMarkerSplitRunPplacer(remove_pct={self.remove_pct})', title=True)
        self.make_output_dirs()

        # Note: These were run on a HPC so pplacer must be called manually (see slurm.sh)
        # They were also split to make it a lot faster.

        with tempfile.TemporaryDirectory() as tmp_dir:
            for i in tqdm(range(self.n_trees_output), total=self.n_trees_output):
                pplacer = Pplacer()

                pplacer_json = os.path.join(self.root_dir, f'pplacer_{i}.json')
                path_tree_tmp = os.path.join(tmp_dir, f'pplacer_{i}.tree')
                pplacer.tog(pplacer_json, path_tree_tmp)
                copy_file(path_tree_tmp, self.output()[str(i)].path, checksum=True)

        return
