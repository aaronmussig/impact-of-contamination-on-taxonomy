import os
import tempfile

import luigi
import pandas as pd
from luigi import LocalTarget
from magna.util.disk import copy_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONTIG_SPLIT
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani_contig_split.c_analyse_fastani import FastAniContigSplitAnalyseFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetTsv
from workflow.util.log import log


class FastAniContigSplitReportResultsFastAni(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    @property
    def root_dir(self):
        return DIR_OUT_FASTANI_CONTIG_SPLIT

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'sp_results': FastAniContigSplitAnalyseFastAni(target_pct=self.target_pct, mash_k=self.mash_k,
                                                           mash_s=self.mash_s)
        }

    def output(self):
        return LocalTargetTsv(
            os.path.join(
                self.root_dir,
                f'results_report_k{self.mash_k}_s{self.mash_s}__p{self.target_pct}.tsv'
            ))

    def run(self):
        log(f'FastAniContigSplitReportResultsFastAni(p={self.target_pct}, mash_k={self.mash_k}, mash_s={self.mash_s})',
            title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading species clustering results')
        df_sp_results = self.input()['sp_results'].maybe_read_cached()
        df_sp_results.set_index('gid', inplace=True)

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()

        log('Parsing results')
        rows = list()
        for gid, row in tqdm(df_max_css.iterrows(), total=len(df_max_css)):

            sp_result_keep = df_sp_results.loc[gid]
            sp_result_disc = df_sp_results.loc[f'{gid}_C']

            expected_tax = df_meta.loc[gid]['gtdb_taxonomy']
            if sp_result_keep.new_sp_rep and isinstance(sp_result_keep.new_sp_rep, str):
                keep_new_tax = df_meta.loc[sp_result_keep.new_sp_rep]['gtdb_taxonomy']
            else:
                keep_new_tax = ''
            if sp_result_disc.new_sp_rep and isinstance(sp_result_disc.new_sp_rep, str):
                disc_new_tax = df_meta.loc[sp_result_disc.new_sp_rep]['gtdb_taxonomy']
            else:
                disc_new_tax = ''

            # Calculate the highest rank that agrees when a new taxonomy is assigned
            if keep_new_tax:
                highest_agree_keep = calc_highest_agree_tax(keep_new_tax, expected_tax)
            else:
                highest_agree_keep = ''
            if disc_new_tax:
                highest_agree_disc = calc_highest_agree_tax(disc_new_tax, expected_tax)
            else:
                highest_agree_disc = ''

            rows.append({
                'gid': gid,
                'expected_tax': expected_tax,
                'keep_tax': keep_new_tax,
                'disc_tax': disc_new_tax,
                'keep_rank_agree': highest_agree_keep,
                'disc_rank_agree': highest_agree_disc,
                'keep_sp_rep': sp_result_keep.new_sp_rep if str(sp_result_keep.new_sp_rep) != 'nan' else '',
                'disc_sp_rep': sp_result_disc.new_sp_rep if str(sp_result_disc.new_sp_rep) != 'nan' else '',
                'keep_ani': sp_result_keep.ani if str(sp_result_keep.ani) != 'nan' else 0.0,
                'disc_ani': sp_result_disc.ani if str(sp_result_disc.ani) != 'nan' else 0.0,
                'keep_af': round(sp_result_keep.af, 4) if str(sp_result_keep.af) != 'nan' else 0.0,
                'disc_af': round(sp_result_disc.af, 4) if str(sp_result_disc.af) != 'nan' else 0.0,
                'keep_type': sp_result_keep.type,
                'disc_type': sp_result_disc.type,
                'keep_same_as_207': sp_result_keep.same,
                'disc_same_as_207': sp_result_disc.same,
                'keep_same_as_disc': str(sp_result_keep.new_sp_rep) == str(sp_result_disc.new_sp_rep)

            })

        df_rows = pd.DataFrame(rows)

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp = os.path.join(tmp_dir, 'out.tsv')
            df_rows.to_csv(path_tmp, sep='\t', index=False)
            log(f'Copying {get_file_size_fmt(path_tmp)} to: {self.output().path}')
            copy_file(path_tmp, self.output().path, checksum=True)
        return


def calc_highest_agree_tax(a, b):
    prev = 'd'
    for cur_a, cur_b in zip(a.split(';'), b.split(';')):
        if cur_a != cur_b:
            return prev
        else:
            prev = cur_a[0]
    return prev
