import os
import subprocess
import tempfile

import luigi
from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_FASTTREE_FULL_TREE_NON_REPS
from workflow.fasttree_full_tree_non_reps.gunc_a_create_batches import FastTreeFullSetNonRepsCreateBatchesGunc
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeFullSetNonRepsCreateBatchesGuncFastTree(LuigiTask):
    congruence = luigi.FloatParameter()
    target_pct = luigi.FloatParameter()


    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, f'c{str(self.congruence)}_p{str(self.target_pct)}')


    def requires(self):
        return {
            'msa': FastTreeFullSetNonRepsCreateBatchesGunc(congruence=self.congruence, target_pct=self.target_pct),
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, 'fasttree_undecorated.tree'))

    def run(self):
        log(f'FastTreeFullSetNonRepsCreateBatchesGuncFastTree(congruence={self.congruence}) (pct={self.target_pct})',
            title=True)
        self.make_output_dirs()

        log('Running FastTree')
        cmd = [
            'FastTreeMP',
            '-wag',
            self.input()['msa'].path,
        ]

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp_out = os.path.join(tmp_dir, 'output.tree')
            path_tmp_stderr = os.path.join(tmp_dir, 'stderr.txt')

            with open(path_tmp_out, 'w') as f_out:
                with open(path_tmp_stderr, 'w') as f_stderr:
                    proc = subprocess.Popen(cmd, stdout=f_out, stderr=f_stderr, encoding='utf-8')
                    proc.communicate()
                    if proc.returncode != 0:
                        raise Exception(f'FastTree failed with return code {proc.returncode}')

            copy_file(path_tmp_out, self.output().path, checksum=True)
            copy_file(path_tmp_stderr, os.path.join(os.path.dirname(self.output().path), 'stderr.txt'), checksum=True)
