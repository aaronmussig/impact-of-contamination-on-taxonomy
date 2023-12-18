import multiprocessing as mp
import os
import tempfile

from luigi import LocalTarget
from magna.util.disk import copy_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, DIR_OUT_MARKER_FAIL
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid_include_mul
from workflow.method.get_msa_from_hits import align_marker
from workflow.model.luigi import LuigiTask
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.collection import iter_batches
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_r207_root


def collect_fna_worker(gid):
    gid_root = get_gid_r207_root(gid)

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}_protein.faa')
    d_faa = dict(read_fasta(path_faa))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # Iterate over each percentage to determine which makers are kept
    markers, hits = get_marker_hits_for_gid_include_mul(d_faa, pfam_th, tigr_th)

    out = list()
    for hit in hits:
        out.append((gid, hit['hit'].gene_id, hit['hit'].hmm_id, hit['n'], hit['seq']))
    return out


class GetFailMarkers(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def load_reps(self):
        df = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()
        return sorted(set(df.index))

    def run(self):
        log('Getting markers for GTDB R207 GUNC failed genomes', title=True)
        self.make_output_dirs()
        os.makedirs(DIR_OUT_MARKER_FAIL, exist_ok=True)

        log('Loading metadata')
        fail_gids = self.load_reps()
        log(f'Loaded: {len(fail_gids):,}')

        batches = list(iter_batches(fail_gids, 1000))
        gids_seen = set()
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir_unaln = os.path.join(tmp_dir, 'unaligned')
            tmp_dir_aln = os.path.join(tmp_dir, 'aligned')

            os.makedirs(tmp_dir_unaln)
            os.makedirs(tmp_dir_aln)

            log(f'Writing to: {tmp_dir_unaln}')

            for batch in tqdm(batches):
                with mp.Pool(processes=mp.cpu_count()) as pool:
                    results = list(pool.imap_unordered(collect_fna_worker, batch))

                    for result in results:
                        for gid, gene_id, hmm_id, n, seq in result:
                            gids_seen.add(gid)
                            path_out = os.path.join(tmp_dir_unaln, f'{hmm_id}.fna')
                            with open(path_out, 'a') as fh:
                                fh.write(f'>{gid}|{gene_id}|{n}\n{seq}\n')

            log(f'Processed: {len(gids_seen):,}')

            log('Aligning markers...')
            gids_seen_2 = set()
            for file_name in os.listdir(tmp_dir_unaln):
                marker = file_name.replace('.fna', '')

                path_aln_out = os.path.join(DIR_OUT_MARKER_FAIL, f'{marker}.faa')

                file_path = os.path.join(tmp_dir_unaln, file_name)
                aligned_markers = align_marker(marker, file_path)

                path_tmp_aligned = os.path.join(tmp_dir_aln, f'{marker}.faa')
                with open(path_tmp_aligned, 'w') as f:
                    for gid, seq in sorted(aligned_markers.items()):
                        gids_seen_2.add(gid[0:15])
                        f.write(f'>{gid}\n{seq}\n')

                if DEBUG:
                    print(path_tmp_aligned)
                else:
                    log(f'Copying {get_file_size_fmt(path_tmp_aligned)} to: {path_aln_out}')
                    copy_file(path_tmp_aligned, path_aln_out, checksum=True)

            log(f'Processed: {len(gids_seen_2):,}')
        if not DEBUG:
            self.write_sentinel()

        return
