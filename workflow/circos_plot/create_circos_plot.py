import multiprocessing as mp
import os
import tempfile

import luigi
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE, MASH_MIN_DIST, DEBUG, DIR_OUT_CIRCOS, R207_BAC120_HMM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.mash_db import GtdbR207MashDb
from workflow.fastani_congruence.a_run_mash_on_failed import FastAniCongruenceRunMash
from workflow.fastani_congruence.b_run_fastani import FastAniCongruenceRunFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.genome_pct_congruence_contigs_removed import GenomePctCongruenceContigsRemoved
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.mash import SketchFile, DistanceFile
from workflow.util.log import log
from workflow.util.paths import get_gid_root

BAC120_MARKERS = frozenset(R207_BAC120_HMM)

class CreateCircosPlot(LuigiTask):

    query_gids = luigi.ListParameter()
    ref_gid = luigi.Parameter()

    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    marker_track = luigi.BoolParameter(default=False)

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def run(self):
        log(f'CreateCircosPlot', title=True)

        out_dir = os.path.join(DIR_OUT_CIRCOS, f'tmp_7')
        log(f'Writing to: {out_dir}')
        os.makedirs(out_dir, exist_ok=True)

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        ref_css = df_max_css.loc[self.ref_gid, 'taxonomic_level']
        ref_source = df_max_css.loc[self.ref_gid, 'source']
        ref_source = GuncRefDb(ref_source)

        log('Loading genomes')
        ref_genome = Genome(self.ref_gid)
        query_fnas = [os.path.join(get_gid_root(gid), f'{gid}.fna') for gid in self.query_gids]

        log('Inferring ordering of contigs (worst to best)')
        df_contig_order = ref_genome.get_gunc_max_css_contig_ranking(ref_css, ref_source)
        d_contig_to_pct_correct = {x.contig: round(x.pct_correct, 2) for x in df_contig_order.itertuples()}
        contig_order = [str(x) for x in df_contig_order['contig']]

        log('Getting marker loci')
        ref_marker_loci = ref_genome.get_marker_loci()


        log('Determining what contigs were put into the C half')
        markers_at_pct_removed, markers_expected, domain_vote, markers_at_pct_removed_c = ref_genome.get_unq_markers_present_at_pct_removed(
            max_css=ref_css, source=ref_source, pct=50.0
        )

        log('Writing re-ordered genome from worst to best contigs')
        path_ref_tmp = os.path.join(out_dir, f'{self.ref_gid}_reordered.fna')
        with open(path_ref_tmp, 'w') as f:
            for contig in contig_order:
                contig_seq = ref_genome.d_fna[contig]
                f.write(f'>{contig}\n{contig_seq}\n')



        log('Writing marker loci file')
        path_labels = os.path.join(out_dir, 'loci.txt')
        pct_values_reported = set()
        with open(path_labels, 'w') as f:
            for marker, lst_loci in ref_marker_loci.items():
                if marker not in BAC120_MARKERS:
                    continue

                for i, (gene_id, locus_from, locus_to) in enumerate(lst_loci):
                    locus_contig = gene_id[:gene_id.rindex('_')]
                    if len(lst_loci) > 1:
                        marker_name = marker + '^{' + str(i + 1) + '/' + str(len(lst_loci)) + '}'
                        # marker_name = '{marker}_{i+1}/{len(lst_loci)}'
                    else:
                        marker_name = marker
                    f.write(f'{locus_contig}\t{locus_from}\t{locus_to}\t{marker_name}\n')

            # # Annotate with congruence value
            # for contig, pct_correct in d_contig_to_pct_correct.items():
            #     if pct_correct not in pct_values_reported:
            #         f.write(f'{contig}\t{0}\t{1}\t{pct_correct}%\tcolor=red\n')
            #         pct_values_reported.add(pct_correct)

        log('Generating command')
        cmd = [
            'mummer2circos',
            '-a', 'nucmer',
            '-l',
            '-lf',
            path_labels,
            '-r',
            path_ref_tmp,
            '-q',
        ]
        cmd.extend(query_fnas)

        path_cmd = os.path.join(out_dir, 'cmd.sh')
        with open(path_cmd, 'w') as f:
            f.write(' '.join(cmd))



        return

        ref_fna = os.path.join(get_gid_root(self.ref_gid), f'{self.ref_gid}.fna')
        query_fnas = [os.path.join(get_gid_root(gid), f'{gid}.fna') for gid in self.query_gids]

        ref_genome = Genome(self.ref_gid)
        ref_marker_loci = ref_genome.get_marker_loci()

        # Get the congruence value of the reference genome


        df_contig_ranking = ref_genome.get_gunc_max_css_contig_ranking(ref_css, ref_source)
        contig_order = list(df_contig_ranking['contig'].values)

        contig_assignments = ref_genome.get_gunc_max_css_contig_assignments(ref_css, ref_source)
        from collections import Counter
        cnts = Counter([contig_assignments[x] for x in contig_order if x in contig_assignments])
        most_common = sorted(cnts.items(), key=lambda x: -x[1])

        # print(f'Reference gid {self.ref_gid} = {df_meta.loc[self.ref_gid, "gtdb_taxonomy"]}')
        # for qry_gid in self.query_gids:
        #     print(f'Query gid {qry_gid} = {df_meta.loc[qry_gid, "gtdb_taxonomy"]}')


        out_dir = os.path.join(DIR_OUT_CIRCOS, f'{self.ref_gid}__' + '_'.join(self.query_gids))
        os.makedirs(out_dir, exist_ok=True)
        log(out_dir)

        # log('Loading ANI')
        # df_ani = self.input()['ani'].maybe_read_cached()
        # df_ani = df_ani[df_ani['query'] == self.ref_gid]
        #
        # log('Loading Mash')
        # df_mash = self.input()['mash'].maybe_read_cached()
        # df_mash = df_mash[df_mash['query'] == self.ref_gid]

        # Save the reference genome re-ordered in contig ranking
        path_ref_tmp = os.path.join(out_dir, f'{self.ref_gid}_reordered.fna')
        with open(path_ref_tmp, 'w') as f:
            for contig in contig_order:
                contig_seq = ref_genome.d_fna[contig]
                f.write(f'>{contig}\n{contig_seq}\n')

        # Write the marker loci file
        path_labels = os.path.join(out_dir, 'loci.txt')
        pct_values_reported = set()
        with open(path_labels, 'w') as f:
            for marker, lst_loci in ref_marker_loci.items():
                # if marker not in BAC120_MARKERS:
                #     continue

                for i, (gene_id, locus_from, locus_to) in enumerate(lst_loci):
                    locus_contig = gene_id[:gene_id.rindex('_')]
                    if len(lst_loci) > 1:
                        marker_name = marker + '^{' + str(i + 1) + '/' + str(len(lst_loci)) + '}'
                        # marker_name = '{marker}_{i+1}/{len(lst_loci)}'
                    else:
                        marker_name = marker
                    f.write(f'{locus_contig}\t{locus_from}\t{locus_to}\t{marker_name}\n')

            # Annotate with congruence value
            for cur_congruence in range(10, 101, 10):
                cur_df_subset = df_contig_ranking[df_contig_ranking['pct_correct'] >= cur_congruence]
                if len(cur_df_subset) == 0:
                    continue
                last_series = df_contig_ranking[df_contig_ranking['pct_correct'] >= cur_congruence].iloc[0]
                pct_value = int(round(float(last_series.pct_correct), 0))
                if pct_value not in pct_values_reported:
                    f.write(f'{last_series.contig}\t{0}\t{1}\t{pct_value}%\tcolor=red\n')
                    pct_values_reported.add(pct_value)

        # Generate the command
        cmd = [
            'mummer2circos',
            '-a', 'promer',
            '-l',
            '-lf',
            path_labels,
            '-r',
            path_ref_tmp,
            '-q',
        ]
        cmd.extend(query_fnas)

        path_cmd = os.path.join(out_dir, 'cmd.sh')
        with open(path_cmd, 'w') as f:
            f.write(' '.join(cmd))

        print(f'Reference gid {self.ref_gid} = {df_meta.loc[self.ref_gid, "gtdb_taxonomy"]}')
        for qry_gid in self.query_gids:
            print(f'Query gid {qry_gid} = {df_meta.loc[qry_gid, "gtdb_taxonomy"]}')

        # print(cmd)

        return

