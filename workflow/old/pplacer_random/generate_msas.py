import multiprocessing as mp
import os
import tempfile

from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_PPLACER_RANDOM_ARC, DIR_OUT_PPLACER_RANDOM_BAC, DIR_OUT_BATCH
from workflow.config import DIR_OUT_SENTINEL
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid
from workflow.method.get_msa_from_hits import get_msa_from_hits
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.pplacer_random.select_contigs_for_pplacer import SelectContigsForPplacer
from workflow.util.collection import iter_batches
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root, get_gid_r207_root
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty

N_CPUS = mp.cpu_count() // 2

PCT = 50

MAX_GIDS = 5000


def worker(job):
    rep, gid_i, gid, pct, contigs_to_remove, domain = job

    # Read the FASTA file
    gid_root = get_gid_r207_root(gid)

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    path_faa = os.path.join(gid_root, 'prodigal',  f'{gid}_protein.faa')

    d_faa = dict(read_fasta(path_faa))

    results_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th, contigs_to_remove)
    try:
        msa = get_msa_from_hits(results_markers, user_domain=domain, marker_dir='/tmp/am/markers' if not DEBUG else None)
    except Exception:
        print(rep, gid_i, gid, pct, domain, contigs_to_remove)
        raise

    return rep, gid_i, gid, pct, msa, domain


class GeneratePplacerRandomMsas(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'bootstrap_contigs': SelectContigsForPplacer(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def get_gid_to_domain(self):
        df = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()
        return df['domain'].to_dict()

    def run(self):
        log('Generating pplacer msa for randomly removed contigs', title=True)
        self.make_output_dirs()

        log('Loading gid to domain')
        d_gid_to_domain = self.get_gid_to_domain()

        log('Loading bootstrap contigs')
        df_contigs = self.input()['bootstrap_contigs'].read() if not DEBUG else self.input()[
            'bootstrap_contigs'].read_cached()

        log('Creating queue')
        queue = list()
        for row in df_contigs.itertuples():

            rep = row.rep
            gid_i = row.gid_i
            gid = row.gid
            pct = row.pct

            domain = d_gid_to_domain[gid]
            queue.append((rep, gid_i, gid, pct, row.contigs_to_remove, domain))

            if DEBUG and len(queue) >= 40:
                break

        log('Writing batch files')
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        batch_paths = list()
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir, exist_ok=True)
            for i, cur_batch in enumerate(iter_batches(queue, n=500 if not DEBUG else 2)):
                cur_batch_path = os.path.join(batch_dir, f'batch_{i}.tsv')
                with open(cur_batch_path, 'w') as f:
                    for rep, gid_i, gid, pct, contigs_to_remove, domain in cur_batch:
                        f.write(f'{rep}\t{gid_i}\t{gid}\t{pct}\t{contigs_to_remove}\t{domain}\n')
                batch_paths.append((cur_batch_path,))

            if not DEBUG:
                log('Submitting batches to RQ...')
                submit_jobs_to_rq(fn=batch_worker, q_args=batch_paths, queue_name=self.fqn)
        else:
            batch_paths = [(os.path.join(batch_dir, x),) for x in os.listdir(batch_dir)]

        if DEBUG:
            log('Starting workers (single-threaded)')
            [batch_worker(x[0]) for x in batch_paths]
        else:
            # Run each of the batches
            log('Waiting until RQ queue is empty...')
            rq_wait_for_queue_empty(self.fqn)
            # self.write_sentinel()

        log('Aggregating results...')
        arc_hits, bac_hits = list(), list()
        results_tld = os.path.join(DIR_OUT_BATCH, 'pplacer_generate_msas_tmp_output')
        for batch_path in batch_paths:
            result_file = os.path.join(results_tld, os.path.basename(batch_path[0]))
            with open(result_file) as f:
                for line in f.readlines():
                    rep, gid_i, gid, pct, msa, domain = line.strip().split('\t')
                    rep = int(rep)
                    gid_i = int(gid_i)
                    pct = int(pct)

                    if domain == 'd__Archaea':
                        arc_hits.append((rep, gid_i, gid, pct, msa))
                    elif domain == 'd__Bacteria':
                        bac_hits.append((rep, gid_i, gid, pct, msa))
                    else:
                        raise ValueError(f'Unexpected domain: {domain}')

        log('Writing results')
        write_hits_to_disk_batched(arc_hits, DIR_OUT_PPLACER_RANDOM_ARC)
        write_hits_to_disk_batched(bac_hits, DIR_OUT_PPLACER_RANDOM_BAC)

        if not DEBUG:
            self.write_sentinel()

        #
        # log(f'Creating MSA for ({len(queue):,})')
        # with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
        #     results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))
        #
        # log('Aggregating into domain-specific MSAs')
        # arc_results, bac_results = list(), list()
        # for rep, gid_i, gid, pct, msa, domain in results:
        #     if domain == 'd__Archaea':
        #         arc_results.append((rep, gid_i, gid, pct, msa))
        #     elif domain == 'd__Bacteria':
        #         bac_results.append((rep, gid_i, gid, pct, msa))
        #     else:
        #         raise ValueError(f'Unknown domain: {domain}')
        #
        # log('Writing to batched files')
        # with tempfile.TemporaryDirectory() as tmp_dir:
        #     for results, domain, target_root in [
        #         (arc_results, 'arc', DIR_OUT_PPLACER_RANDOM_ARC),
        #         (bac_results, 'bac', DIR_OUT_PPLACER_RANDOM_BAC)
        #     ]:
        #         for i, batch in enumerate(iter_batches(results, MAX_GIDS)):
        #             path = os.path.join(tmp_dir, f'{domain}_pplacer_{i}.msa')
        #             with open(path, 'w') as fh:
        #                 for rep, gid_i, gid, pct, msa in batch:
        #                     fh.write(f'>{rep}__{gid_i}__{gid}__{pct}\n{msa}\n')
        #             target_dir = os.path.join(target_root, str(i))
        #             os.makedirs(target_dir, exist_ok=True)
        #             path_out = os.path.join(target_dir, f'{domain}_pplacer.msa')
        #             copy_file(path, path_out, checksum=True)
        #
        # if not DEBUG:
        #     self.write_sentinel()


def batch_worker(batch_path):
    log(f'Reading batch file: {batch_path}')
    queue = list()
    with open(batch_path) as f:
        for line in f.readlines():
            rep, gid_i, gid, pct, contigs_to_remove, domain = line.strip().split('\t')
            queue.append((int(rep), int(gid_i), gid, int(pct), frozenset(contigs_to_remove.split('|')), domain))
    log(f'Found {len(queue):,} entries in batch file')

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

    tld = os.path.join(DIR_OUT_BATCH, 'pplacer_generate_msas_tmp_output')
    os.makedirs(tld, exist_ok=True)
    path_out = os.path.join(DIR_OUT_BATCH, 'pplacer_generate_msas_tmp_output', os.path.basename(batch_path))

    with open(path_out, 'w') as f:
        for rep, gid_i, gid, pct, msa, domain in results:
            f.write(f'{rep}\t{gid_i}\t{gid}\t{pct}\t{msa}\t{domain}\n')
    return


def write_hits_to_disk_batched(hits, root):
    log(f'Writing {len(hits):,} hits to disk')
    for i, batch in enumerate(iter_batches(hits, MAX_GIDS)):
        out_root = os.path.join(root, str(i))
        os.makedirs(out_root, exist_ok=True)
        path = os.path.join(out_root, 'input.fasta')
        with open(path, 'w') as fh:
            for rep, gid_i, gid, pct, msa in batch:
                fh.write(f'>{rep}__{gid_i}__{gid}__{pct}\n{msa}\n')
    return
