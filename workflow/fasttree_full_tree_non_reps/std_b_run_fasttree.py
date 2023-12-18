import os
import subprocess
import tempfile

from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_FASTTREE_FULL_TREE_NON_REPS
from workflow.fasttree_full_tree_non_reps.std_a_create_batches import FastTreeFullSetNonRepsCreateBatchesStandard
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeFullSetNonRepsCreateBatchesStandardFastTree(LuigiTask):

    def requires(self):
        return {
            'msa': FastTreeFullSetNonRepsCreateBatchesStandard(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, 'standard', 'fasttree_undecorated.tree'))

    def run(self):
        log(f'FastTreeFullSetNonRepsCreateBatchesGuncFastTree()',
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
            log(f'stderr: {path_tmp_stderr}')
            log(f'output: {path_tmp_out}')

            with open(path_tmp_out, 'w') as f_out:
                with open(path_tmp_stderr, 'w') as f_stderr:
                    proc = subprocess.Popen(cmd, stdout=f_out, stderr=f_stderr, encoding='utf-8')
                    proc.communicate()
                    if proc.returncode != 0:
                        raise Exception(f'FastTree failed with return code {proc.returncode}')

            copy_file(path_tmp_out, self.output().path, checksum=True)
            copy_file(path_tmp_stderr, os.path.join(os.path.dirname(self.output().path), 'stderr.txt'), checksum=True)
