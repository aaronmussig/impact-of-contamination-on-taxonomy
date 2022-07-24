import os
import tempfile

from magna.util.disk import md5sum, move_file
from magna.util.web import download_file

from workflow.config import DIR_OUT_EXTERNAL
from workflow.model.luigi import LocalTargetTree, LuigiTask
from workflow.util.log import log


class GtdbTreeArcR207(LuigiTask):
    """External GTDB Archaeal tree, available for download from the GTDB website."""

    def output(self):
        return LocalTargetTree(os.path.join(DIR_OUT_EXTERNAL, 'gtdb.r207.ar122.tree'))

    def run(self):
        log('Loading GTDB R207 Archaeal tree...', title=True)
        self.make_output_dirs()

        path_ar122 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_r207.tree'
        ar122_md5 = 'e6e493882a067e18d1d6fd50be89c047'

        with tempfile.TemporaryDirectory() as tmpdir:
            path_ar122_tmp = os.path.join(tmpdir, 'tree.tree')
            download_file(path_ar122, path_ar122_tmp)

            if md5sum(path_ar122_tmp) != ar122_md5:
                raise ValueError('Checksum mismatch for ar122 tree')

            move_file(path_ar122_tmp, self.output().path, checksum=True)


class GtdbTreeBacR207(LuigiTask):
    """External GTDB Bacterial tree, available for download from the GTDB website."""

    def output(self):
        return LocalTargetTree(os.path.join(DIR_OUT_EXTERNAL, 'gtdb.r207.bac120.tree'))

    def run(self):
        log('Loading GTDB R207 Bacterial tree...', title=True)
        self.make_output_dirs()

        path_bac120 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_r207.tree'
        bac120_md5 = '12fab8ba201d200aada9c8d61dafda4b'

        with tempfile.TemporaryDirectory() as tmpdir:
            path_bac120_tmp = os.path.join(tmpdir, 'tree.tree')
            download_file(path_bac120, path_bac120_tmp)

            if md5sum(path_bac120_tmp) != bac120_md5:
                raise ValueError('Checksum mismatch for bac120 tree')

            move_file(path_bac120_tmp, self.output().path, checksum=True)
