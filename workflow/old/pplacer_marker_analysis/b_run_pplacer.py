import multiprocessing as mp
import os
import tempfile

import luigi
from magna.util.disk import copy_file

from workflow.config import PPLACER_R207_ARC_PATH, \
    PPLACER_R207_BAC_PATH, DIR_OUT_PPLACER_MARKER_ANALYSIS
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.model.pplacer import Pplacer
from workflow.pplacer_marker_analysis.a_create_batches import PplacerMarkerAnalysisCreateBatches
from workflow.util.log import log


class PplacerMarkerAnalysisRunPplacer(LuigiTask):
    # Required for both
    domain = luigi.ChoiceParameter(choices=['arc', 'bac'])
    target_pct = luigi.FloatParameter()

    # Required for GUNC-informed
    congruence = luigi.FloatParameter(default=-1)

    # Required for random
    random = luigi.BoolParameter(default=False)
    batch_id = luigi.IntParameter(default=-1)

    def requires(self):
        return {
            'msa': PplacerMarkerAnalysisCreateBatches(
                domain=self.domain,
                target_pct=self.target_pct,
                congruence=self.congruence,
                random=self.random,
                batch_id=self.batch_id
            )
        }

    def output(self):
        if self.random:
            target_dir = os.path.join(
                DIR_OUT_PPLACER_MARKER_ANALYSIS,
                f'random_{self.domain}_p{self.target_pct}_b{self.batch_id}'
            )
        else:
            target_dir = os.path.join(
                DIR_OUT_PPLACER_MARKER_ANALYSIS,
                f'gunc_{self.domain}_p{self.target_pct}_c{self.congruence}'
            )
        return LocalTargetTree(os.path.join(target_dir, 'pplacer.tree'))

    def run(self):
        if self.random:
            log(f'RANDOM: PplacerMarkerAnalysisRunPplacer('
                f'pct={self.target_pct}, batch_id={self.batch_id}, '
                f'domain={self.domain})', title=True)
        else:
            log(f'GUNC: PplacerMarkerAnalysisRunPplacer('
                f'congruence={self.congruence}, pct={self.target_pct}, '
                f'domain={self.domain})', title=True)
        self.make_output_dirs()

        out_dir = os.path.dirname(self.output().path)
        msa_path = self.input()['msa']['msa'].path

        if self.domain == 'arc':
            path_refpkg = PPLACER_R207_ARC_PATH
        else:
            path_refpkg = PPLACER_R207_BAC_PATH

        pplacer = Pplacer()
        pplacer_json_out = os.path.join(out_dir, 'pplacer.json')
        pplacer_out = os.path.join(out_dir, 'pplacer.out')

        log(f'Writing pplacer log to: {pplacer_out}')
        log(f'Writing pplacer json to: {pplacer_json_out}')

        pplacer.run(min(mp.cpu_count(), 64), 'wag', path_refpkg, pplacer_json_out, msa_path, pplacer_out)

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tree_tmp = os.path.join(tmp_dir, 'pplacer.tree')
            pplacer.tog(pplacer_json_out, path_tree_tmp)
            copy_file(path_tree_tmp, self.output().path, checksum=True)

        return
