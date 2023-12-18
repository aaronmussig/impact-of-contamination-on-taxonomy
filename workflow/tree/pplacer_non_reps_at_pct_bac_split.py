import multiprocessing as mp
import os
import tempfile
from collections import defaultdict
from typing import Tuple, List

from luigi import LocalTarget
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, PPLACER_R207_ARC_PATH, \
    DIR_OUT_PPLACER_ARC_NON_REPS, DIR_OUT_PPLACER_BAC_NON_REPS, PPLACER_R207_BAC_PATH, DIR_OUT_SENTINEL
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_msa_for_failed_at_pct import GetMsaForFailedAtPct
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.model.pplacer import Pplacer
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait

OUTPUT_DIRECTORY = DIR_OUT_PPLACER_BAC_NON_REPS


def msa_collect_worker(gid) -> List[Tuple[str, int, str]]:
    out = list()

    gid_root = get_gid_root(gid)
    path_msa = os.path.join(gid_root, 'msa_at_cutoff_values.tsv')
    last_msa = None

    with open(path_msa) as f:
        f.readline()
        for line in f.readlines():
            result_gid, result_pct, result_msa = line.strip().split('\t')

            if set(result_msa) == {'-'}:
                print(f'Warning: No information for {result_gid} at {result_pct}')
                continue

            if result_pct == '0':
                last_msa = result_msa
            else:
                if last_msa != result_msa:
                    out.append((result_gid, int(result_pct), result_msa))
                    last_msa = result_msa

    return out


def collect_all_msas(gids):
    if DEBUG:
        gids = gids[:10]

    unq_gids = set()
    unq_msas = set()

    out = defaultdict(list)

    if DEBUG:
        for gid in gids:
            for result_gid, result_pct, result_msa in msa_collect_worker(gid):
                out[result_pct].append((result_gid, result_msa))
    else:
        with mp.Pool(processes=min(mp.cpu_count(), 40)) as pool:
            for results in tqdm(pool.imap_unordered(msa_collect_worker, gids), total=len(gids)):
                for result_gid, result_pct, result_msa in results:
                    out[result_pct].append((result_gid, result_msa))
                    unq_gids.add(result_gid)
                    unq_msas.add(result_msa)

    print(f'Found: {len(unq_gids):,} unique gids and {len(unq_msas):,} unique msas')
    return out

def pplacer_worker(output_dir):
    log(f'Running pplacer: {output_dir}')

    path_msa = os.path.join(output_dir, 'msa.fasta')

    pplacer = Pplacer()
    pplacer_json_out = os.path.join(output_dir, 'pplacer.json')
    pplacer_out = os.path.join(output_dir, 'pplacer.out')

    pplacer.run(min(mp.cpu_count(), 64), 'wag', PPLACER_R207_BAC_PATH, pplacer_json_out, path_msa, pplacer_out)

    path_tree_file = os.path.join(output_dir, 'output.tree')

    pplacer.tog(pplacer_json_out, path_tree_file)

    return


class PplacerNonRepsAtPctBacSplit(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            '_msa_at_pct': GetMsaForFailedAtPct(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Running pplacer for bac (split)', title=True)

        self.make_output_dirs()
        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        log('Subsetting to non-rep bact genomes')
        df_merged = df_merged[(df_merged['domain'] == 'd__Bacteria')]
        log(f'Found {len(df_merged):,} bacterial genomes')

        queue = list()

        log('Getting all MSAs that differ from the baseline')
        with tempfile.TemporaryDirectory(prefix='gunc-chim_') as tmp_dir:
            # Write all MSAs to a temp directory
            d_pct_to_msa = collect_all_msas(df_merged.index)

            # Write to files
            for pct, result_list in d_pct_to_msa.items():
                path_tmp = os.path.join(tmp_dir, f'msa_{pct}.fasta')
                with open(path_tmp, 'w') as f:
                    for gid, msa in sorted(result_list):
                        f.write(f'>{gid}_{pct}\n{msa}\n')

                # Copy to output directory
                cur_pct_dir = os.path.join(OUTPUT_DIRECTORY, str(pct))
                os.makedirs(cur_pct_dir, exist_ok=True)
                path_msa_out = os.path.join(cur_pct_dir, 'msa.fasta')
                copy_file(path_tmp, path_msa_out, checksum=True)
                queue.append((cur_pct_dir, ))

        # Run each of the batches
        log('Submitting batches to RQ...')
        rq_and_wait(job_id=self.fqn, fn=pplacer_worker, q_args=queue, queue_name=self.fqn)

        return
