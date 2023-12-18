import multiprocessing as mp
import os
import tempfile

import luigi
from magna.util.disk import copy_file

from workflow.config import DIR_OUT_EXTERNAL, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.model.luigi import LuigiTask, BaseLocalTarget
from workflow.model.mash import SketchFile
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class GtdbR207MashDb(LuigiTask):
    """Mash database of all GTDB R207 genomes."""

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return BaseLocalTarget(os.path.join(
            DIR_OUT_EXTERNAL,
            f'gtdb_r207_mash_db_k{self.mash_k}_s{self.mash_s}.msh'
        ))

    def run(self):
        log(f'Creating GTDB R207 Mash reference database (k={self.mash_k}, s={self.mash_s})', title=True)
        self.make_output_dirs()

        df_meta = self.input()['meta'].maybe_read_cached()
        rep_gids = set(df_meta[df_meta['gtdb_representative'] == 't'].index)
        log(f'Found {len(rep_gids):,} GTDB representative genomes')

        d_gids = dict()
        for gid in sorted(rep_gids):
            d_gids[gid] = os.path.join(get_gid_root(gid), f'{gid}.fna')
            if DEBUG and len(d_gids) > 10:
                break

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_sketch_tmp = os.path.join(tmp_dir, 'sketch.msh')
            log(f'Temporary file path: {path_sketch_tmp}')
            SketchFile(genomes=d_gids, path=path_sketch_tmp, cpus=mp.cpu_count(), k=self.mash_k, s=self.mash_s)
            log(f'Copying to: {self.output().path} ')
            copy_file(path_sketch_tmp, self.output().path, checksum=True)

        return




class GtdbR207MashDbAll(LuigiTask):
    """Mash database of all GTDB R207 genomes."""

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return BaseLocalTarget(os.path.join(
            DIR_OUT_EXTERNAL,
            f'gtdb_r207_mash_db_all_k{self.mash_k}_s{self.mash_s}.msh'
        ))

    def run(self):
        log(f'Creating GTDB R207 Mash reference database (all) (k={self.mash_k}, s={self.mash_s})', title=True)
        self.make_output_dirs()

        df_meta = self.input()['meta'].maybe_read_cached()
        all_gids = set(df_meta.index)
        log(f'Found {len(all_gids):,} GTDB genomes')

        d_gids = dict()
        for gid in sorted(all_gids):
            d_gids[gid] = os.path.join(get_gid_root(gid), f'{gid}.fna')
            if DEBUG and len(d_gids) > 10:
                break

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_sketch_tmp = os.path.join(tmp_dir, 'sketch.msh')
            log(f'Temporary file path: {path_sketch_tmp}')
            SketchFile(genomes=d_gids, path=path_sketch_tmp, cpus=mp.cpu_count(), k=self.mash_k, s=self.mash_s)
            log(f'Copying to: {self.output().path} ')
            copy_file(path_sketch_tmp, self.output().path, checksum=True)

        return
