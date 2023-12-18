import os
import tempfile

from luigi import LocalTarget
from magna.util.disk import md5sum, move_file
from magna.util.web import download_file

from workflow.config import DIR_OUT_EXTERNAL, PATH_R207_AR53_MASK, PATH_R207_BAC120_MASK
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class GtdbMaskArcR207(LuigiTask):
    """External GTDB Archaeal tree, available for download from the GTDB website."""

    def output(self):
        return LocalTarget(PATH_R207_AR53_MASK)

    def run(self):
        log('Loading GTDB R207 Archaeal mask...', title=True)
        self.make_output_dirs()

        path_ar122 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/ar53_msa_mask_r207.txt'
        ar122_md5 = '18952c7861dc8f0db7d0b24892495e91'

        with tempfile.TemporaryDirectory() as tmpdir:
            path_ar122_tmp = os.path.join(tmpdir, 'tree.tree')
            download_file(path_ar122, path_ar122_tmp)

            if md5sum(path_ar122_tmp) != ar122_md5:
                raise ValueError('Checksum mismatch for ar122 mask')

            move_file(path_ar122_tmp, self.output().path, checksum=True)


class GtdbMaskBacR207(LuigiTask):
    """External GTDB Bacterial tree, available for download from the GTDB website."""

    def output(self):
        return LocalTarget(PATH_R207_BAC120_MASK)

    def run(self):
        log('Loading GTDB R207 Bacterial mask...', title=True)
        self.make_output_dirs()

        path_bac120 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/bac120_msa_mask_r207.txt'
        bac120_md5 = '38f050479ac768453bc7c60d106ecca3'

        with tempfile.TemporaryDirectory() as tmpdir:
            path_bac120_tmp = os.path.join(tmpdir, 'tree.tree')
            download_file(path_bac120, path_bac120_tmp)

            if md5sum(path_bac120_tmp) != bac120_md5:
                raise ValueError('Checksum mismatch for bac120 mask')

            move_file(path_bac120_tmp, self.output().path, checksum=True)
