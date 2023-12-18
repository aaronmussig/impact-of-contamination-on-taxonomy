import os
import subprocess
import tempfile

import luigi
from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_MARKER_TRIMMING_GUNC
from workflow.marker_trimming.gunc_a_create_msa import MarkerTrimmingGuncCreateMsa
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class MarkerTrimmingGuncRunFastTree(LuigiTask):
    target_pct = luigi.FloatParameter()
    replicate_id = luigi.IntParameter()

    @property
    def root_dir(self):
        if self.replicate_id == 0:
            return os.path.join(DIR_OUT_MARKER_TRIMMING_GUNC, f'p{str(self.target_pct)}')
        else:
            return os.path.join(DIR_OUT_MARKER_TRIMMING_GUNC, f'p{str(self.target_pct)}', 'fasttree_replicates')

    def requires(self):
        return {
            'msa': MarkerTrimmingGuncCreateMsa(target_pct=self.target_pct),
        }

    def output(self):
        return LocalTarget(os.path.join(self.root_dir, f'fasttree_undecorated_r{self.replicate_id}.tree'))

    def run(self):
        log(f'MarkerTrimmingGuncRunFastTree(pct={self.target_pct}, replicate={self.replicate_id})', title=True)
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

            path_stderr_out = os.path.join(self.root_dir, f'fasttree_stderr_r{self.replicate_id}.txt')

            copy_file(path_tmp_out, self.output().path, checksum=True)
            copy_file(path_tmp_stderr, path_stderr_out, checksum=True)
