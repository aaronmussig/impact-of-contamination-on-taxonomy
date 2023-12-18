import multiprocessing as mp
import os
from collections import defaultdict

import pandas as pd
from chiton.fastani import fastani
from tqdm import tqdm

from workflow.config import DIR_OUT_V2_FASTANI_INTERSP, DIR_OUT_BATCH, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait

BATCH_SIZE = 2500

GIDS_NEW_SP = {'GCA_001509115.1', 'GCA_001695755.1', 'GCA_002728285.1', 'GCA_002731855.1', 'GCA_008668585.1',
               'GCA_008668795.1', 'GCA_009493725.1', 'GCA_011523145.1', 'GCA_016707075.1', 'GCA_017394825.1',
               'GCA_017465765.1', 'GCA_017515185.1', 'GCA_018056875.1', 'GCA_018239885.1', 'GCA_018363345.1',
               'GCA_900759525.1', 'GCA_900761055.1', 'GCA_900765305.1', 'GCA_900765645.1', 'GCA_902528895.1',
               'GCA_902593295.1', 'GCA_903846615.1', 'GCA_903931905.1', 'GCA_905200745.1', 'GCA_905214645.1',
               'GCF_000698005.1', 'GCF_002929465.1'}
D_GIDS_CHG_SP = {'GCA_002291775.1': 'GCA_002712565.1',
                 'GCA_002703565.1': 'GCA_002690085.1',
                 'GCA_003487585.1': 'GCA_002712565.1',
                 'GCA_007096555.1': 'GCA_007095785.1',
                 'GCA_011523095.1': 'GCA_002690995.1',
                 'GCA_013213925.1': 'GCF_001550135.1',
                 'GCA_016939315.1': 'GCF_000785515.1',
                 'GCA_018370815.1': 'GCF_000478885.1',
                 'GCA_018383955.1': 'GCA_900554835.1',
                 'GCA_900548685.1': 'GCF_019041915.1',
                 'GCA_900549575.1': 'GCF_003459245.1',
                 'GCA_900756095.1': 'GCA_900542085.1',
                 'GCA_900760075.1': 'GCF_001405375.1',
                 'GCA_902363945.1': 'GCA_900543275.1',
                 'GCA_902388545.1': 'GCA_900543275.1',
                 'GCA_902539685.1': 'GCA_902619825.1',
                 'GCA_902568005.1': 'GCA_902619825.1',
                 'GCA_902573475.1': 'GCA_902546975.1',
                 'GCA_902619155.1': 'GCA_902555695.1',
                 'GCA_905193315.1': 'GCA_900546775.1',
                 'GCA_905208535.1': 'GCA_900544115.1',
                 'GCA_905212495.1': 'GCA_900541065.1',
                 'GCA_905214135.1': 'GCA_902399975.1',
                 'GCA_905214255.1': 'GCF_003436275.1',
                 'GCA_905215365.1': 'GCF_003459245.1',
                 'GCF_000964075.1': 'GCF_000813765.1',
                 'GCF_001812365.1': 'GCF_000311745.1',
                 'GCF_001815585.1': 'GCF_000311745.1',
                 'GCF_002035805.1': 'GCF_008692785.1',
                 'GCF_002835625.1': 'GCF_900475885.1',
                 'GCF_002912425.1': 'GCF_003048675.2',
                 'GCF_002913635.1': 'GCF_003048675.2',
                 'GCF_003462845.1': 'GCA_900543275.1',
                 'GCF_003495765.1': 'GCF_001457655.1',
                 'GCF_008269775.1': 'GCF_000307025.1',
                 'GCF_009493905.1': 'GCF_015074785.1',
                 'GCF_009494015.1': 'GCF_000157935.1',
                 'GCF_009649595.1': 'GCF_000016825.1',
                 'GCF_014284905.1': 'GCF_000742135.1',
                 'GCF_015667745.1': 'GCF_001405375.1',
                 'GCF_016806285.1': 'GCF_900478295.1',
                 'GCF_016806405.1': 'GCF_900478295.1',
                 'GCF_017151305.1': 'GCF_002208095.1',
                 'GCF_018499505.1': 'GCF_001073155.1',
                 'GCF_902845755.1': 'GCF_900478295.1'}


class V2FastAniInterSpeciesGetValuesForFailedToGenus(LuigiTask):
    """
    This script finds the pairwise ANI between all species representatives in
    each genus that contains at least one failed genome (bacteria only).
    """

    @property
    def root_dir(self):
        return DIR_OUT_V2_FASTANI_INTERSP

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, f'failed_to_genus.h5'))

    def run(self):
        log(f'V2FastAniInterSpeciesGetValuesForFailedToGenus', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading MaxCSS')
        df_css = self.input()['css'].maybe_read_cached()
        failed_gids = frozenset(df_css.index)
        failed_bac_gids = failed_gids.intersection(df_meta[df_meta['domain'] == 'd__Bacteria'].index)
        failed_bac_gids_no_change = frozenset(failed_bac_gids - set(D_GIDS_CHG_SP.keys()) - GIDS_NEW_SP)

        log('Finding the genera for each failed genome')
        failed_genera = frozenset(df_meta[df_meta.index.isin(failed_gids)]['genus'].unique())
        bac_genera = frozenset(df_meta[df_meta['domain'] == 'd__Bacteria']['genus'].unique())
        bac_failed_genera = bac_genera.intersection(failed_genera)
        log(f'Found {len(bac_failed_genera):,} bacterial genera with failed genomes')

        log('Finding the species representatives for each genus')
        sp_reps_for_failed_genera = frozenset(
            df_meta[(df_meta['genus'].isin(bac_failed_genera)) & (df_meta['gtdb_representative'] == 't')].index)
        log(f'Found {len(sp_reps_for_failed_genera):,} species representatives for failed genera')

        d_gid_to_sp = df_meta[df_meta['gtdb_representative'] == 't']['species'].to_dict()
        d_sp_to_gid = {v: k for k, v in d_gid_to_sp.items()}

        log('Generating a mapping of genus to species in each')
        d_genus_to_sp_reps = defaultdict(set)
        for gid, row in df_meta[df_meta.index.isin(sp_reps_for_failed_genera)].iterrows():
            d_genus_to_sp_reps[row['genus']].add(gid)

        log('Splitting into three cases')

        # 1. For those that changed species, compare the ANI between the two species reps that changed
        sp_change_comparisons = dict()
        for fail_gid, new_gid in D_GIDS_CHG_SP.items():
            cur_fail_sp = df_meta.loc[fail_gid, 'species']
            cur_new_sp = df_meta.loc[new_gid, 'species']
            sp_change_comparisons[(cur_fail_sp, cur_new_sp)] = (d_sp_to_gid[cur_fail_sp], d_sp_to_gid[cur_new_sp])

        # 2. For those that formed a new species cluster, compare the ANI between all species in the genera
        sp_new_cluster_comparisons = dict()
        for fail_gid in GIDS_NEW_SP:
            cur_fail_genus = df_meta.loc[fail_gid, 'genus']
            cur_fail_genus_sp_reps = d_genus_to_sp_reps[cur_fail_genus]
            sp_new_cluster_comparisons[cur_fail_genus] = cur_fail_genus_sp_reps

        # 3. For those that did not chage, compare ANI between all speices in the genera
        sp_no_change_comparisons = dict()
        for fail_gid in failed_bac_gids_no_change:
            cur_fail_genus = df_meta.loc[fail_gid, 'genus']
            cur_fail_genus_sp_reps = d_genus_to_sp_reps[cur_fail_genus]
            sp_no_change_comparisons[cur_fail_genus] = cur_fail_genus_sp_reps

        # Check if the batches have already been created
        batch_dir = os.path.join(DIR_OUT_BATCH, self.fqn)
        if not os.path.isdir(batch_dir):
            log(f'Creating batch files: {batch_dir}')
            os.makedirs(batch_dir, exist_ok=True)

            log('Creating array of pairwise comparisons to perform')
            comparisons = create_compairson_array(d_genus_to_sp_reps)

            log('Splitting comparisons into batches')
            batches = list()
            for batch_i, batch in enumerate(iter_batches(comparisons, BATCH_SIZE)):
                batch_path = os.path.join(batch_dir, f'batch_{batch_i}.tsv')
                with open(batch_path, 'w') as f:
                    for x, y in batch:
                        f.write(f'{x}\t{y}\n')
                batches.append((batch_path,))

            log('Submitting jobs to RQ')
            rq_and_wait(job_id=self.fqn, fn=batch_worker, q_args=batches, queue_name=self.fqn)

        log('Loading batch directory')
        batch_paths = load_from_batch_dir(batch_dir)
        log(f'Found {len(batch_paths):,} batch files')

        # if True:
        #     [batch_worker(x) for x in batch_paths]

        # if True:
        #     log('Submitting jobs to RQ')
        #     rq_and_wait(job_id=self.fqn, fn=batch_worker, q_args=[(x,) for x in batch_paths], queue_name=self.fqn)

        log('Concatenating results')
        dfs = list()
        for path in batch_paths:
            # if len(dfs) > 5 and DEBUG:
            #     break
            try:
                dfs.append(pd.read_csv(path.replace('.tsv', '.result.tsv'), sep='\t'))
            except Exception:
                continue

        df_raw = pd.concat(dfs, ignore_index=True)

        log('Scoring results')



        log('Processing results into the genus')
        df_genus = convert_df_raw(df_raw, d_genus_to_sp_reps, sp_change_comparisons, sp_new_cluster_comparisons, sp_no_change_comparisons)

        return


def convert_df_raw(df_raw, d_genus_to_sp_reps, sp_change_comparisons, sp_new_cluster_comparisons, sp_no_change_comparisons):
    # Index the rows for faster access
    d_q_v_r = dict()
    for row in df_raw.itertuples():
        d_q_v_r[(row.query, row.ref)] = (float(row.ani), float(row.af))
    df_raw = None

    # Case 1.
    case_1_ani, case_1_af = list(), list()
    for (_, _), (gid1, gid2) in sp_change_comparisons.items():
        try:
            ani, af = d_q_v_r.get((gid1, gid2), d_q_v_r[(gid2, gid1)])
            if ani > 0 and af > 0.5:  # always true
                case_1_ani.append(ani)
                case_1_af.append(af)
        except Exception:
            continue

    # Case 2.
    case_2_ani, case_2_af = list(), list()
    for genus, sp_reps in sp_new_cluster_comparisons.items():
        comparisons = create_compairson_array(d_genus_to_sp_reps, genus)
        for gid1, gid2 in comparisons:
            try:
                ani, af = d_q_v_r.get((gid1, gid2), d_q_v_r[(gid2, gid1)])
                if ani > 0 and af > 0.5:
                    case_2_ani.append(ani)
                    case_2_af.append(af)
            except Exception:
                continue

    # Case 3.
    case_3_ani, case_3_af = list(), list()
    for genus, sp_reps in sp_no_change_comparisons.items():
        comparisons = create_compairson_array(d_genus_to_sp_reps, genus)
        for gid1, gid2 in comparisons:
            try:
                ani, af = d_q_v_r.get((gid1, gid2), d_q_v_r[(gid2, gid1)])
                if ani > 0 and af > 0.5:
                    case_3_ani.append(ani)
                    case_3_af.append(af)
            except Exception:
                continue

    return


def mp_worker(job):
    qry_gid, ref_gid = job
    rows = run_fastani_on_jobs(ref_gid, qry_gid)
    return rows


def batch_worker(path):
    queue = list()
    with open(path) as f:
        for line in f.readlines():
            x, y = line.strip().split('\t')
            queue.append((x, y))

    log(f'Processing batch: {path}')
    all_rows = list()
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap_unordered(mp_worker, queue), total=len(queue)))
        [all_rows.extend(x) for x in results]

    df = pd.DataFrame(all_rows)
    path_out = path.replace('.tsv', '.result.tsv')
    df.to_csv(path_out, index=False, sep='\t')
    return


def load_from_batch_dir(batch_dir):
    out = list()
    batch_paths = [os.path.join(batch_dir, x) for x in os.listdir(batch_dir) if x.endswith('tsv')]
    for batch_path in batch_paths:
        out.append(batch_path)
    return out


def create_compairson_array(d_genus_to_sp_reps, target_genus=None):
    out = list()
    for genus, sp_reps in d_genus_to_sp_reps.items():
        if target_genus is not None and target_genus != genus:
            continue
        sp_reps_lst = list(sp_reps)
        cur_comparisons = set()
        for i in range(len(sp_reps_lst)):
            for j in range(i + 1, len(sp_reps_lst)):
                cur_comparisons.add((sp_reps_lst[i], sp_reps_lst[j]))
        out.extend(sorted(cur_comparisons))

    return out


def run_fastani_on_jobs(ref_gid: str, qry_gid: str):
    ref_gid_paths = {ref_gid: os.path.join(get_gid_root(ref_gid), f'{ref_gid}.fna')}
    qry_gid_paths = {qry_gid: os.path.join(get_gid_root(qry_gid), f'{qry_gid}.fna')}

    out = list()

    ani = fastani(query=list(qry_gid_paths.values()),
                  reference=list(ref_gid_paths.values()),
                  cpus=1,
                  single_execution=False,
                  bidirectional=True,
                  show_progress=False)

    # Prepare the ANI results for output
    d_ani = ani.as_dict()
    for qry_gid in qry_gid_paths.keys():
        q_key = qry_gid_paths[qry_gid]
        for ref_gid in ref_gid_paths.keys():
            r_key = ref_gid_paths[ref_gid]
            qvr = d_ani[q_key][r_key]
            rvq = d_ani[r_key][q_key]

            if qvr is not None and rvq is not None:
                ani = max(qvr.ani, rvq.ani)
                af = max(qvr.align_frac, rvq.align_frac)
            elif qvr is not None and rvq is None:
                ani = qvr.ani
                af = qvr.align_frac
            elif qvr is None and rvq is not None:
                ani = rvq.ani
                af = rvq.align_frac
            else:
                ani = 0
                af = 0

            out.append({
                'query': qry_gid,
                'ref': ref_gid,
                'ani': ani,
                'af': af,
            })
    return out
