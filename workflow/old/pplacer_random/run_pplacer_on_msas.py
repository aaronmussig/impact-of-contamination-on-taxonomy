import multiprocessing as mp
import os

from magna.util.disk import move_file

from workflow.config import DIR_OUT_PPLACER_RANDOM_ARC, DIR_OUT_PPLACER_RANDOM_BAC, PPLACER_R207_ARC_PATH, \
    PPLACER_R207_BAC_PATH, DEBUG
from workflow.config import DIR_OUT_SENTINEL
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.pplacer import Pplacer
from workflow.pplacer_random.generate_msas import GeneratePplacerRandomMsas
from workflow.util.log import log
from workflow.util.rq import submit_jobs_to_rq, rq_wait_for_queue_empty


def enqueue_from_dir(dir_name):
    out = list()
    for batch in os.listdir(dir_name):
        path_msa = os.path.join(dir_name, batch, 'input.fasta')
        # new_path = os.path.join(dir_name, batch, 'input.fasta')
        # move_file(path_msa, new_path, checksum=True)
        out.append(path_msa)
    log(f'Enqueued {len(out):,} msas from {dir_name}')
    return out


def run_pplacer(output_directory, path_msa, path_refpkg):
    pplacer = Pplacer()
    pplacer_json_out = os.path.join(output_directory, 'pplacer.json')
    pplacer_out = os.path.join(output_directory, 'pplacer.out')

    pplacer.run(min(mp.cpu_count(), 64), 'wag', path_refpkg, pplacer_json_out, path_msa, pplacer_out)

    path_tree_file = os.path.join(output_directory, 'pplacer_random.tree')
    pplacer.tog(pplacer_json_out, path_tree_file)
    return


class RunPplacerRandomOnMsas(LuigiTask):

    def requires(self):
        return {
            '_pplacer_msas': GeneratePplacerRandomMsas(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running pplacer on MSAs (random)', title=True)
        self.make_output_dirs()

        redis_queue = list()

        arc_input = enqueue_from_dir(DIR_OUT_PPLACER_RANDOM_ARC)
        for path in arc_input:
            out_directory = os.path.dirname(path)
            redis_queue.append((out_directory, path, PPLACER_R207_ARC_PATH))

        bac_input = enqueue_from_dir(DIR_OUT_PPLACER_RANDOM_BAC)
        for path in bac_input:
            out_directory = os.path.dirname(path)
            redis_queue.append((out_directory, path, PPLACER_R207_BAC_PATH))

        log('Submitting batches to RQ...')
        submit_jobs_to_rq(fn=run_pplacer, q_args=redis_queue, queue_name=self.fqn)

        log('Waiting until RQ queue is empty...')
        rq_wait_for_queue_empty(self.fqn)

        if not DEBUG:
            self.write_sentinel()
