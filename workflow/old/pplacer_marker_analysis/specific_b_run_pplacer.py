import multiprocessing as mp
import os
import tempfile

import luigi
from magna.util.disk import copy_file

from workflow.config import PPLACER_R207_ARC_PATH, \
    PPLACER_R207_BAC_PATH, DIR_OUT_PPLACER_MARKER_ANALYSIS
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.model.pplacer import Pplacer
from workflow.pplacer_marker_analysis.specific_a_create_msa import PplacerMarkerAnalysisSpecificCreateMsa
from workflow.util.log import log


class PplacerMarkerAnalysisSpecificRunPplacer(LuigiTask):
    pplacer_domain = luigi.ChoiceParameter(choices=['arc', 'bac'])
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()
    gids = luigi.ListParameter()

    def requires(self):
        return {
            'msa': PplacerMarkerAnalysisSpecificCreateMsa(
                pplacer_domain=self.pplacer_domain,
                target_pct=self.target_pct,
                congruence=self.congruence,
                gids=self.gids
            )
        }

    @property
    def root_dir(self) -> str:
        return os.path.join(
            DIR_OUT_PPLACER_MARKER_ANALYSIS,
            f'specific_{self.pplacer_domain}_p{self.target_pct}_c{self.congruence}__{"_".join(sorted(self.gids))}'
        )

    def output(self):
        return LocalTargetTree(os.path.join(self.root_dir, 'pplacer.tree'))

    def run(self):
        log(f'PplacerMarkerAnalysisSpecificRunPplacer('
            f'congruence={self.congruence}, pct={self.target_pct}, '
            f'domain={self.pplacer_domain}, gids={self.gids})', title=True)
        self.make_output_dirs()

        msa_path = self.input()['msa']['msa'].path

        if self.pplacer_domain == 'arc':
            path_refpkg = PPLACER_R207_ARC_PATH
        else:
            path_refpkg = PPLACER_R207_BAC_PATH

        pplacer = Pplacer()
        pplacer_json_out = os.path.join(self.root_dir, 'pplacer.json')
        pplacer_out = os.path.join(self.root_dir, 'pplacer.out')

        log(f'Writing pplacer log to: {pplacer_out}')
        log(f'Writing pplacer json to: {pplacer_json_out}')

        pplacer.run(min(mp.cpu_count(), 64), 'wag', path_refpkg, pplacer_json_out, msa_path, pplacer_out)

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tree_tmp = os.path.join(tmp_dir, 'pplacer.tree')
            pplacer.tog(pplacer_json_out, path_tree_tmp)
            copy_file(path_tree_tmp, self.output().path, checksum=True)

        return
