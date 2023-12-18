import os
import subprocess
import tempfile

import luigi
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_FASTTREE_RANDOMLY_REMOVE_MARKERS
from workflow.fasttree_randomly_remove_markers.a_create_msa import FastTreeRandomlyRemoveMarkersCreateMsa
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.util.log import log

"""
This was run on a HPC instead
"""


class FastTreeRandomlyRemoveMarkersRunFastTree(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)
    replicate_id = luigi.IntParameter()

    @property
    def root_dir(self):
        return os.path.join(
            DIR_OUT_FASTTREE_RANDOMLY_REMOVE_MARKERS,
            f'pct_{self.target_pct}',
            f'replicate_{self.replicate_id}'
        )

    def requires(self):
        return {
            'msa': FastTreeRandomlyRemoveMarkersCreateMsa(target_pct=self.target_pct, replicate_id=self.replicate_id),
        }

    def output(self):
        return LocalTargetTree(os.path.join(self.root_dir, f'fasttree_undecorated_unrooted.tree'))

    def run(self):
        log(f'FastTreeRandomlyRemoveMarkersRunFastTree(pct={self.target_pct})', title=True)
        self.make_output_dirs()

        log('Running FastTree')
        cmd = [
            'FastTreeMP',
            '-wag',
            self.input()['msa']['msa'].path,
        ]

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp_out = os.path.join(tmp_dir, 'output.tree')
            path_tmp_stderr = os.path.join(tmp_dir, 'stderr.txt')
            log(f'Writing FastTree output to {path_tmp_out}')
            log(f'Writing FastTree stderr to {path_tmp_stderr}')

            with open(path_tmp_out, 'w') as f_out:
                with open(path_tmp_stderr, 'w') as f_stderr:
                    proc = subprocess.Popen(cmd, stdout=f_out, stderr=f_stderr, encoding='utf-8')
                    proc.communicate()
                    if proc.returncode != 0:
                        raise Exception(f'FastTree failed with return code {proc.returncode}')

            path_stderr_out = os.path.join(self.root_dir, f'fasttree_undecorated_unrooted.stderr')

            copy_file(path_tmp_out, self.output().path, checksum=True)
            copy_file(path_tmp_stderr, path_stderr_out, checksum=True)
