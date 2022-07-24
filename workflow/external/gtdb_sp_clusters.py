import os
import tempfile

import pandas as pd
from gunc.gunc_database import download_file
from magna.util.disk import md5sum

from workflow.config import DIR_OUT_EXTERNAL
from workflow.model.luigi import LocalTargetHdf5, LuigiTask
from workflow.util.log import log


class GtdbSpClustersR207(LuigiTask):
    """External GTDB species clusters file, available for download from the GTDB website."""

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_EXTERNAL, 'gtdb.r207.sp_clusters.h5'))

    def run(self):
        log('Creating GTDB R207 species clusters file...', title=True)
        self.make_output_dirs()

        path = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/sp_clusters_r207.tsv'
        md5 = 'a8952b2090e7d092acc49a6528bf85fc'

        # Download into a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            log('Downloading file...')
            path_tmp = os.path.join(tmpdir, 'sp_clusters_r207.tsv')
            download_file(path, path_tmp)

            log('Verifying checksums...')
            if md5sum(path_tmp) != md5:
                raise ValueError('Checksum mismatch for species clusters file')

            log('Loading dataframe...')
            df = pd.read_csv(path_tmp, sep='\t', index_col=False, low_memory=False)

        # Modify table
        df.rename(columns={
            'Representative genome': 'rep_genome',
            'GTDB species': 'species',
            'GTDB taxonomy': 'taxonomy',
            'ANI circumscription radius': 'ani_radius',
            'Mean intra-species ANI': 'ani_mean',
            'Min intra-species ANI': 'ani_min',
            'Mean intra-species AF': 'af_mean',
            'Min intra-species AF': 'af_min',
            'No. clustered genomes': 'n_genomes',
            'Clustered genomes': 'clustered_genomes'
        }, inplace=True)
        df.drop(columns=['clustered_genomes'], inplace=True)

        # Save to disk
        self.save_hdf(df, index='species')
