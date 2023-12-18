import multiprocessing as mp
import os
import tempfile
from typing import Tuple, List

from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, PPLACER_R207_ARC_PATH, \
    DIR_OUT_PPLACER_ARC_NON_REPS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_msa_for_failed_at_pct import GetMsaForFailedAtPct
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.model.pplacer import Pplacer
from workflow.util.log import log
from workflow.util.paths import get_gid_root

OUTPUT_DIRECTORY = DIR_OUT_PPLACER_ARC_NON_REPS


def msa_collect_worker(gid) -> List[Tuple[str, int, str]]:
    out = list()

    gid_root = get_gid_root(gid)
    path_msa = os.path.join(gid_root, 'msa_at_cutoff_values.tsv')
    last_msa = None

    with open(path_msa) as f:
        f.readline()
        for line in f.readlines():
            result_gid, result_pct, result_msa = line.strip().split('\t')

            if result_pct == '0':
                last_msa = result_msa
            else:
                if last_msa != result_msa:
                    out.append((result_gid, int(result_pct), result_msa))
                    last_msa = result_msa

    return out


def collect_all_msas(gids, path_msa):
    if DEBUG:
        gids = gids[:10]

    unq_gids = set()
    unq_msas = set()

    with open(path_msa, 'w') as f:
        if DEBUG:
            for gid in gids:
                for result_gid, result_pct, result_msa in msa_collect_worker(gid):
                    f.write(f'>{result_gid}_{result_pct}\n{result_msa}\n')
        else:
            with mp.Pool(processes=min(mp.cpu_count(), 40)) as pool:
                for results in tqdm(pool.imap_unordered(msa_collect_worker, gids), total=len(gids)):
                    for result_gid, result_pct, result_msa in results:
                        f.write(f'>{result_gid}_{result_pct}\n{result_msa}\n')
                        unq_gids.add(result_gid)
                        unq_msas.add(result_msa)

    print(f'Found: {len(unq_gids):,} unique gids and {len(unq_msas):,} unique msas')
    return


class PplacerNonRepsAtPctArc(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            '_msa_at_pct': GetMsaForFailedAtPct(),
        }

    def output(self):
        return LocalTargetTree(os.path.join(OUTPUT_DIRECTORY, 'output.tree'))

    def run(self):
        log('Running pplacer for non arc representatives at pct values', title=True)
        self.make_output_dirs()

        self.make_output_dirs()
        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        log('Subsetting to non-rep archaeal genomes')
        df_merged = df_merged[(df_merged['domain'] == 'd__Archaea')]
        log(f'Found {len(df_merged):,} archaeal genomes')

        path_msa_out = os.path.join(OUTPUT_DIRECTORY, 'msa.fasta')

        log('Getting all MSAs that differ from the baseline')
        with tempfile.TemporaryDirectory(prefix='gunc-chim_') as tmp_dir:
            # Write all MSAs to a temp directory
            user_msa_file = os.path.join(tmp_dir, 'msa.fasta')
            collect_all_msas(df_merged.index, user_msa_file)
            copy_file(user_msa_file, path_msa_out, checksum=True)

        log('Running pplacer')
        pplacer = Pplacer()
        pplacer_json_out = os.path.join(OUTPUT_DIRECTORY, 'pplacer.json')
        pplacer_out = os.path.join(OUTPUT_DIRECTORY, 'pplacer.out')

        pplacer.run(min(mp.cpu_count(), 64), 'wag', PPLACER_R207_ARC_PATH, pplacer_json_out, path_msa_out, pplacer_out)

        path_tree_file = self.output().path

        pplacer.tog(pplacer_json_out, path_tree_file)

        return
