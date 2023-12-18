import multiprocessing as mp
import os
import tempfile
from pathlib import Path

import luigi
from gtdblib.tree.bootstrap_merge import bootstrap_merge_replicates
from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_MARKER_TRIMMING_GUNC, DEBUG
from workflow.marker_trimming.gunc_c_root_fasttree import MarkerTrimmingGuncRootFastTree
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class MarkerTrimmingGuncCalcBootstrapFastTree(LuigiTask):
    target_pct = luigi.FloatParameter()

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_MARKER_TRIMMING_GUNC, f'p{str(self.target_pct)}')

    def requires(self):
        d_requires = {
            'tree': MarkerTrimmingGuncRootFastTree(target_pct=self.target_pct, replicate_id=0),
        }
        for i in range(1, 101):
            d_requires[f'tree_{i}'] = MarkerTrimmingGuncRootFastTree(target_pct=self.target_pct, replicate_id=i)
        return d_requires

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'fasttree_undecorated_rooted_support.tree'))

    def run(self):
        log(f'MarkerTrimmingGuncCalcBootstrapFastTree(pct={self.target_pct}', title=True)
        self.make_output_dirs()

        path_input = Path(self.input()['tree'].path)
        rep_paths = list()
        for i in range(1, 101):
            key = f'tree_{i}'
            if key in self.input():
                rep_paths.append(Path(self.input()[key].path))

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            path_out_tmp = tmp_dir / 'output.tree'

            bootstrap_merge_replicates(path_input, path_out_tmp, rep_paths, mp.cpu_count())

            if not DEBUG:
                log(f'Copying file to: {self.output().path}')
                copy_file(str(path_out_tmp), self.output().path, checksum=True)

        return
