import multiprocessing as mp
import os
import tempfile
from pathlib import Path

from gtdblib.tree.bootstrap_merge import bootstrap_merge_replicates
from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DEBUG, DIR_OUT_FASTTREE_FULL_TREE_AND_FAILED_NO_TRIM
from workflow.fasttree_full_tree_and_failed_no_trim.c_root_fasttree import FastTreeFullTreeAndFailedNoTrimRoot
from workflow.model.luigi import LuigiTask
from workflow.util.log import log

import sys
sys.setrecursionlimit(100000000)


class FastTreeFullTreeAndFailedNoTrimBootstrapRoot(LuigiTask):

    @property
    def root_dir(self):
        return DIR_OUT_FASTTREE_FULL_TREE_AND_FAILED_NO_TRIM

    def requires(self):
        d_requires = {
            'tree': FastTreeFullTreeAndFailedNoTrimRoot(replicate_id=0),
        }
        for i in range(1, 101):
            d_requires[f'tree_{i}'] = FastTreeFullTreeAndFailedNoTrimRoot(replicate_id=i)
        return d_requires

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'fasttree_undecorated_rooted_support.tree'))

    def run(self):
        log(f'FastTreeFullTreeAndFailedNoTrimBootstrapRoot', title=True)
        self.make_output_dirs()

        path_input = Path(self.input()['tree'].path)
        rep_paths = list()
        for i in range(1, 101):
            rep_paths.append(Path(self.input()[f'tree_{i}'].path))

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            path_out_tmp = tmp_dir / 'output.tree'

            bootstrap_merge_replicates(path_input, path_out_tmp, rep_paths,mp.cpu_count())

            if not DEBUG:
                log(f'Copying file to: {self.output().path}')
                copy_file(str(path_out_tmp), self.output().path, checksum=True)

        return
