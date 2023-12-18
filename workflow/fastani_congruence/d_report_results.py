import os
import tempfile

import luigi
import pandas as pd
from luigi import LocalTarget
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani_congruence.b_run_fastani import FastAniCongruenceRunFastAni
from workflow.fastani_congruence.c_analyse_fastani import FastAniCongruenceAnalyseFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.genome_pct_congruence_contigs_removed import GenomePctCongruenceContigsRemoved
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.luigi import LuigiTask
from workflow.util.log import log
from collections import Counter



class FastAniCongruenceReportResults(LuigiTask):
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            # 'ani': FastAniCongruenceRunFastAni(congruence=self.congruence, target_pct=self.target_pct,
            #                                    mash_k=self.mash_k, mash_s=self.mash_s),
            'meta': GtdbMetadataR207(),
            # 'sp_clusters': GtdbSpClustersR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'contigs_removed': GenomePctCongruenceContigsRemoved(congruence=self.congruence, target_pct=self.target_pct),
            'sp_results': FastAniCongruenceAnalyseFastAni(congruence=self.congruence, target_pct=self.target_pct, mash_k=self.mash_k, mash_s=self.mash_s)
        }

    def output(self):
        return LocalTarget(
            os.path.join(
                DIR_OUT_FASTANI_CONGRUENCE,
                f'results_report_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}.tsv'
            ))

    def run(self):
        log(f'FastAniCongruenceReportResults(p={self.target_pct}, c={self.congruence}, mash_k={self.mash_k}, mash_s={self.mash_s})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        bac_gids = frozenset(df_meta[df_meta['domain'] == 'd__Bacteria'].index)

        log('Loading species clustering results')
        df_sp_results = self.input()['sp_results'].maybe_read_cached()

        df_sp_results_wrong = df_sp_results[df_sp_results['same'] == False]
        df_sp_results_wrong = df_sp_results_wrong.copy()

        log('Determining % of genome removed')
        df_contigs_removed = self.input()['contigs_removed'].maybe_read_cached()

        pct_removed_vals = list()

        original_tax = list()
        sp_rep_tax = list()
        for row in df_sp_results_wrong.itertuples():
            contigs_removed_row = df_contigs_removed.loc[row.gid]
            pct_removed = float(contigs_removed_row['pct_removed'])
            # contigs_removed = frozenset(contigs_removed_row['contigs_removed'].split('|'))
            pct_removed_vals.append(pct_removed)

            gid_tax = df_meta.loc[row.gid]['gtdb_taxonomy']
            if row.new_sp_rep and isinstance(row.new_sp_rep, str) :
                new_sp_rep_tax = df_meta.loc[row.new_sp_rep]['gtdb_taxonomy']
            else:
                new_sp_rep_tax = ''
            original_tax.append(gid_tax)
            sp_rep_tax.append(new_sp_rep_tax)

        df_sp_results_wrong['pct_removed'] = pct_removed_vals
        df_sp_results_wrong['original_tax'] = original_tax
        df_sp_results_wrong['sp_rep_tax'] = sp_rep_tax

        # sp_rep no_ani no_af
        cnt_types = Counter(df_sp_results_wrong['type'].values)

        n_total = len(df_sp_results)
        n_wrong = len(df_sp_results_wrong)
        n_correct = n_total - n_wrong
        print(f'Number correct: {n_correct:,}/{n_total:,} = {n_correct/n_total:.2%}')
        print(f'Number wrong: {n_wrong:,}/{n_total:,} = {n_wrong/n_total:.2%}')
        [print(f'{k} = {v} / {n_total} = {v / n_total:.2%}') for k, v in cnt_types.items()]

        df_sp_results_wrong.to_csv('/tmp/fastani_c50_p50.tsv', sep='\t')
        return


