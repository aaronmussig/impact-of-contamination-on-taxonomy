import os
import tempfile

import pandas as pd
from magna.util.disk import untar, md5sum
from magna.util.web import download_file

from workflow.config import DIR_OUT_EXTERNAL
from workflow.model.luigi import LocalTargetHdf5
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class GtdbMetadataR207(LuigiTask):
    """External GTDB metadata file, available for download from the GTDB website."""

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_EXTERNAL, 'gtdb_r207_metadata.h5'))

    def run(self):
        log('Creating GTDB R207 metadata file...', title=True)
        self.make_output_dirs()

        path_ar53 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz'
        ar53_md5 = '610b08c1c1da4e3539a6c88bfa81be5e'

        path_bac120 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz'
        bac120_md5 = 'fb1911b5aa12c098fd25387c4296b157'

        # Download into a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            log('Downloading files...')
            path_ar53_tmp = os.path.join(tmpdir, 'metadata_ar53.tar.gz')
            path_bac120_tmp = os.path.join(tmpdir, 'metadata_bac120.tar.gz')
            download_file(path_ar53, path_ar53_tmp)
            download_file(path_bac120, path_bac120_tmp)

            log('Extracting files...')
            untar(path_ar53_tmp, tmpdir)
            untar(path_bac120_tmp, tmpdir)

            log('Verifying checksums...')
            if md5sum(path_ar53_tmp) != ar53_md5:
                raise ValueError('Checksum mismatch for ar53 metadata file')
            if md5sum(path_bac120_tmp) != bac120_md5:
                raise ValueError('Checksum mismatch for bac120 metadata file')

            log('Loading metadata...')
            df_ar53 = pd.read_csv(os.path.join(tmpdir, 'ar53_metadata_r207.tsv'), sep='\t', index_col=False,
                                  low_memory=False)
            df_bac120 = pd.read_csv(os.path.join(tmpdir, 'bac120_metadata_r207.tsv'), sep='\t', index_col=False,
                                    low_memory=False)

        log('Merging dataframes')
        df = pd.concat([df_ar53, df_bac120], ignore_index=True)
        df['accession'] = df['accession'].apply(lambda x: x[3:])
        df['domain'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[0])
        df['phylum'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[1])
        df['class'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[2])
        df['order'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[3])
        df['family'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[4])
        df['genus'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[5])
        df['species'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[6])
        df.rename(columns={"accession": "gid"}, inplace=True)

        # df['ncbi_molecule_count'] = df['ncbi_molecule_count'].astype(str)
        coerce_cols = [
            'ncbi_molecule_count',
            'ncbi_spanned_gaps',
            'ncbi_species_taxid',
            'ncbi_taxid',
            'ncbi_total_gap_length',
            'ncbi_total_length',
            'ncbi_ungapped_length',
            'ncbi_unspanned_gaps',
        ]
        for col in coerce_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        keep_cols = {'gid', 'gtdb_taxonomy', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species',
                     'ambiguous_bases', 'checkm_completeness', 'checkm_contamination', 'checkm_marker_count',
                     'checkm_marker_lineage', 'checkm_marker_set_count', 'checkm_strain_heterogeneity', 'coding_bases',
                     'coding_density', 'contig_count', 'gc_count', 'gc_percentage', 'genome_size',
                     'gtdb_genome_representative', 'gtdb_representative', 'gtdb_taxonomy', 'gtdb_type_designation',
                     'gtdb_type_designation_sources', 'gtdb_type_species_of_genus',
                     'ncbi_taxonomy', 'ncbi_taxonomy_unfiltered'}
        for col in df.columns:
            if col not in keep_cols:
                df.drop(col, axis=1, inplace=True)

        self.save_hdf(df, index='gid')


class GtdbMetadataR207Full(LuigiTask):
    """External GTDB metadata file, available for download from the GTDB website."""

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_EXTERNAL, 'gtdb_r207_metadata_full.h5'))

    def run(self):
        log('Creating GTDB R207 metadata file...', title=True)
        self.make_output_dirs()

        path_ar53 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz'
        ar53_md5 = '610b08c1c1da4e3539a6c88bfa81be5e'

        path_bac120 = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz'
        bac120_md5 = 'fb1911b5aa12c098fd25387c4296b157'

        # Download into a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            log('Downloading files...')
            path_ar53_tmp = os.path.join(tmpdir, 'metadata_ar53.tar.gz')
            path_bac120_tmp = os.path.join(tmpdir, 'metadata_bac120.tar.gz')
            download_file(path_ar53, path_ar53_tmp)
            download_file(path_bac120, path_bac120_tmp)

            log('Extracting files...')
            untar(path_ar53_tmp, tmpdir)
            untar(path_bac120_tmp, tmpdir)

            log('Verifying checksums...')
            if md5sum(path_ar53_tmp) != ar53_md5:
                raise ValueError('Checksum mismatch for ar53 metadata file')
            if md5sum(path_bac120_tmp) != bac120_md5:
                raise ValueError('Checksum mismatch for bac120 metadata file')

            log('Loading metadata...')
            df_ar53 = pd.read_csv(os.path.join(tmpdir, 'ar53_metadata_r207.tsv'), sep='\t', index_col=False,
                                  low_memory=False)
            df_bac120 = pd.read_csv(os.path.join(tmpdir, 'bac120_metadata_r207.tsv'), sep='\t', index_col=False,
                                    low_memory=False)

        log('Merging dataframes')
        df = pd.concat([df_ar53, df_bac120], ignore_index=True)
        df['accession'] = df['accession'].apply(lambda x: x[3:])
        df['domain'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[0])
        df['phylum'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[1])
        df['class'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[2])
        df['order'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[3])
        df['family'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[4])
        df['genus'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[5])
        df['species'] = df['gtdb_taxonomy'].apply(lambda x: x.split(';')[6])
        df.rename(columns={"accession": "gid"}, inplace=True)

        # df['ncbi_molecule_count'] = df['ncbi_molecule_count'].astype(str)
        coerce_cols = [
            'ncbi_molecule_count',
            'ncbi_spanned_gaps',
            'ncbi_species_taxid',
            'ncbi_taxid',
            'ncbi_total_gap_length',
            'ncbi_total_length',
            'ncbi_ungapped_length',
            'ncbi_unspanned_gaps',
        ]
        for col in coerce_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        self.save_hdf(df, index='gid')
