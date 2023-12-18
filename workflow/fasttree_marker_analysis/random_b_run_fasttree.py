import os
import subprocess
import tempfile

import luigi
from luigi import LocalTarget
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM
from workflow.fasttree_marker_analysis.random_a_create_batches import FastTreeMarkerAnalysisBatchesRandom
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class FastTreeMarkerAnalysisRunRandomFastTree(LuigiTask):
    target_pct = luigi.FloatParameter()
    batch_id = luigi.IntParameter()

    def requires(self):
        return {
            'msa': FastTreeMarkerAnalysisBatchesRandom(target_pct=self.target_pct, batch_id=self.batch_id),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, f'{self.batch_id}_{self.target_pct}',
                                        'fasttree_undecorated.tree'))

    def run(self):
        log(f'Running FastTree bac (random) (batch={self.batch_id}) (pct={self.target_pct})', title=True)
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
            log(f'Writing progress (temp) to: {path_tmp_stderr}')

            with open(path_tmp_out, 'w') as f_out:
                with open(path_tmp_stderr, 'w') as f_stderr:
                    proc = subprocess.Popen(cmd, stdout=f_out, stderr=f_stderr, encoding='utf-8')
                    proc.communicate()
                    if proc.returncode != 0:
                        raise Exception(f'FastTree failed with return code {proc.returncode}')

            copy_file(path_tmp_out, self.output().path, checksum=True)
            copy_file(path_tmp_stderr, os.path.join(os.path.dirname(self.output().path), 'stderr.txt'), checksum=True)
