import os
import tempfile

from magna.util.disk import md5sum, move_file, untar
from magna.util.web import download_file
from tqdm import tqdm

from workflow.config import PATH_R207_AR53_MSA, R207_AR53_HMM, PATH_R207_BAC120_MSA, R207_BAC120_HMM, \
    PATH_R207_AR53_REP_MSA, PATH_R207_BAC120_REP_MSA
from workflow.external.gtdb_mask import GtdbMaskArcR207, GtdbMaskBacR207
from workflow.model.luigi import LuigiTask, LocalTargetFasta
from workflow.util.fasta import read_fasta
from workflow.util.log import log


class GtdbMsaArcNonRepsR207(LuigiTask):
    """Re-create the GTDB R207 MSA for all archaeal genomes"""

    def output(self):
        return LocalTargetFasta(PATH_R207_AR53_MSA)

    def requires(self):
        return {
            'mask': GtdbMaskArcR207(),
        }

    def run(self):
        log('Creating GTDB MSA for all R207 archaeal genomes', title=True)
        self.make_output_dirs()

        log('Reading mask')
        with open(self.input()['mask'].path) as f:
            mask = [x == '1' for x in f.read().strip()]

        path_remote = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_all/ar53_msa_marker_genes_all_r207.tar.gz'
        remote_md5 = 'e28b12457ca57a47139e258cc42f518a'

        d_gid_to_msa = dict()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_archive = os.path.join(tmpdir, 'archive.tar.gz')

            log('Downloading archive')
            download_file(path_remote, tmp_archive)

            log('Validating archive')
            if md5sum(tmp_archive) != remote_md5:
                raise ValueError('Checksum mismatch for file')

            log('Extracting archive')
            untar(tmp_archive, tmpdir)

            log('Extracting markers')
            for marker in tqdm(R207_AR53_HMM):
                path_marker = os.path.join(tmpdir, f'ar53_r207_all_{marker}.faa')

                for gid, seq in read_fasta(path_marker).items():
                    if gid not in d_gid_to_msa:
                        d_gid_to_msa[gid] = list()
                    d_gid_to_msa[gid].append(seq)

            # Concatenate the MSA and write to file
            log(f'Creating MSA for {len(d_gid_to_msa):,} genomes')
            path_msa_tmp = os.path.join(tmpdir, 'msa_output.faa')
            with open(path_msa_tmp, 'w') as f:
                for gid, lst_seqs in tqdm(sorted(d_gid_to_msa.items())):
                    seq = ''.join(lst_seqs)
                    masked_seq = ''.join([x for x, m in zip(seq, mask) if m])
                    f.write(f'>{gid}\n{masked_seq}\n')

            # Move
            log('Moving file')
            move_file(path_msa_tmp, self.output().path, checksum=True)


class GtdbMsaArcR207(LuigiTask):
    """Re-create the GTDB R207 MSA for all archaeal genomes"""

    def output(self):
        return LocalTargetFasta(PATH_R207_AR53_REP_MSA)

    def run(self):
        log('Creating GTDB MSA for all R207 archaeal rep genomes', title=True)
        self.make_output_dirs()

        path_remote = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/ar53_msa_reps_r207.tar.gz'
        remote_md5 = 'f49e8d6b3984b0f05d77209f78d09be0'

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_archive = os.path.join(tmpdir, 'archive.tar.gz')

            log('Downloading archive')
            download_file(path_remote, tmp_archive)

            log('Validating archive')
            if md5sum(tmp_archive) != remote_md5:
                raise ValueError('Checksum mismatch for file')

            log('Extracting archive')
            untar(tmp_archive, tmpdir)
            path_msa_tmp = os.path.join(tmpdir, 'ar53_msa_reps_r207.faa')

            # Move
            log('Moving file')
            move_file(path_msa_tmp, self.output().path, checksum=True)


class GtdbMsaBacNonRepsR207(LuigiTask):
    """Re-create the GTDB R207 MSA for all bac genomes"""

    def output(self):
        return LocalTargetFasta(PATH_R207_BAC120_MSA)

    def requires(self):
        return {
            'mask': GtdbMaskBacR207(),
        }

    def run(self):
        log('Creating GTDB MSA for all R207 bacterial genomes', title=True)
        self.make_output_dirs()

        log('Reading mask')
        with open(self.input()['mask'].path) as f:
            mask = [x == '1' for x in f.read().strip()]

        path_remote = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_all/bac120_msa_marker_genes_all_r207.tar.gz'
        remote_md5 = '97def7dcf7b3f6ed6c1a834416a9c587'

        d_gid_to_msa = dict()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_archive = os.path.join(tmpdir, 'archive.tar.gz')

            log('Downloading archive')
            download_file(path_remote, tmp_archive)

            log('Validating archive')
            if md5sum(tmp_archive) != remote_md5:
                raise ValueError('Checksum mismatch for file')

            log('Extracting archive')
            untar(tmp_archive, tmpdir)

            log('Extracting markers')
            for marker in tqdm(R207_BAC120_HMM):
                path_marker = os.path.join(tmpdir, f'bac120_r207_all_{marker}.faa')

                for gid, seq in read_fasta(path_marker).items():
                    if gid not in d_gid_to_msa:
                        d_gid_to_msa[gid] = list()
                    d_gid_to_msa[gid].append(seq)

            # Concatenate the MSA and write to file
            log(f'Creating MSA for {len(d_gid_to_msa):,} genomes')
            path_msa_tmp = os.path.join(tmpdir, 'msa_output.faa')
            with open(path_msa_tmp, 'w') as f:
                for gid, lst_seqs in tqdm(sorted(d_gid_to_msa.items())):
                    seq = ''.join(lst_seqs)
                    masked_seq = ''.join([x for x, m in zip(seq, mask) if m])
                    f.write(f'>{gid}\n{masked_seq}\n')

            # Move
            log('Moving file')
            move_file(path_msa_tmp, self.output().path, checksum=True)


class GtdbMsaBacR207(LuigiTask):
    """Re-create the GTDB R207 MSA for all bac genomes"""

    def output(self):
        return LocalTargetFasta(PATH_R207_BAC120_REP_MSA)

    def run(self):
        log('Creating GTDB MSA for all R207 bacterial rep genomes', title=True)
        self.make_output_dirs()

        path_remote = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/bac120_msa_reps_r207.tar.gz'
        remote_md5 = '2b90624e7d666198836af396359e8ee4'

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_archive = os.path.join(tmpdir, 'archive.tar.gz')

            log('Downloading archive')
            download_file(path_remote, tmp_archive)

            log('Validating archive')
            if md5sum(tmp_archive) != remote_md5:
                raise ValueError('Checksum mismatch for file')

            log('Extracting archive')
            untar(tmp_archive, tmpdir)
            path_msa_tmp = os.path.join(tmpdir, 'bac120_msa_reps_r207.faa')

            # Move
            log('Moving file')
            move_file(path_msa_tmp, self.output().path, checksum=True)
