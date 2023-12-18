import multiprocessing as mp
import os
import tempfile

from Bio import SeqIO
from luigi import LocalTarget
from magna.gunc import read_contig_assignments_tsv
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_SENTINEL, R207_MARKERS, DIR_PFAM_33, DIR_TIGRFAM_15
from workflow.external.gtdb_mask import GtdbMaskArcR207, GtdbMaskBacR207
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.contig_removal import get_taxonomy_by_majority_vote_gunc, contigs_to_remove_from_gunc
from workflow.method.get_genome_domain_from_markers import get_genome_domain_from_markers
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid
from workflow.method.get_msa_from_hits import get_msa_from_hits
from workflow.model.luigi import LuigiTask
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root


def cache_markers_to_dir(cache_dir):
    log(f'Caching to: {cache_dir}')
    for marker in tqdm(R207_MARKERS):
        if marker.startswith('PF'):
            path = os.path.join(DIR_PFAM_33, 'individual_hmms', f'{marker}.hmm')
        elif marker.startswith('TIGR'):
            path = os.path.join(DIR_TIGRFAM_15, 'individual_hmms', f'{marker}.HMM')
        else:
            raise ValueError(f'Unknown marker: {marker}')

        path_target = os.path.join(cache_dir, f'{marker}.hmm')
        if not os.path.isfile(path_target):
            copy_file(path, path_target, checksum=True)

    return


class GetMsaForFailedAtPct(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            '_mask_arc': GtdbMaskArcR207(),
            '_mask_bac': GtdbMaskBacR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def run(self):
        log('Generating MSA for those gids that failed', title=True)

        self.make_output_dirs()
        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        log('Caching pfam/tigrfam to local disk')
        with tempfile.TemporaryDirectory(prefix='gunc-chim') as tmp_dir:

            log('Caching markers to disk')
            cache_markers_to_dir(tmp_dir)

            queue = list()
            for gid, row in df_merged.iterrows():
                queue.append([
                    gid,
                    row['source'],
                    row['domain'],
                    row['taxonomic_level'],
                    tmp_dir
                ])

                if DEBUG and len(queue) > 3:
                    break

            if DEBUG:
                results = [run_on_gid(x) for x in queue]
            else:
                with mp.Pool(processes=mp.cpu_count()) as pool:
                    results = list(tqdm(pool.imap_unordered(run_on_gid, queue), total=len(queue)))

        log('Saving dataframe')
        if not DEBUG:
            self.write_sentinel()


def run_on_gid(job):
    gid, gunc_source, domain, max_css, marker_dir = job

    # Read the FASTA file
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')

    path_out = os.path.join(gid_root, 'msa_at_cutoff_values.tsv')

    if os.path.isfile(path_out):
        return

    with open(gid_fna) as f:
        d_fna = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    if gunc_source == 'gtdb':
        path_contig_assign = os.path.join(gid_root, 'gunc_r95/gunc_output', f'{gid}.contig_assignments.tsv')
        gunc_domain = domain
    elif gunc_source == 'progenomes':
        path_contig_assign = os.path.join(gid_root, 'gunc_pro/gunc_output', f'{gid}.contig_assignments.tsv')
        if domain == 'd__Bacteria':
            gunc_domain = '2 Bacteria'
        elif domain == 'd__Archaea':
            gunc_domain = '2157 Archaea'
        else:
            raise ValueError(f'Unknown domain: {domain}')
    else:
        raise Exception(f'Unknown gunc source: {gunc_source}')
    df_contig_assign = read_contig_assignments_tsv(path_contig_assign)

    # Determine the percentage values at which this genome can have contigs removed
    taxon, tax_level = get_taxonomy_by_majority_vote_gunc(df_contig_assign, max_css, gunc_domain)
    d_pct_to_contigs_to_remove = contigs_to_remove_from_gunc(d_fna, df_contig_assign, taxon, tax_level)

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_faa = dict(read_fasta(path_faa))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # Iterate over each percentage to determine which makers are kept
    with tempfile.TemporaryDirectory() as tmp_dir:
        path_tmp = os.path.join(tmp_dir, 'out.tsv')
        with open(path_tmp, 'w') as f:
            f.write('gid\tpct\tmsa\n')

            base_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th)
            base_msa = get_msa_from_hits(base_markers, marker_dir=marker_dir)
            base_domain = get_genome_domain_from_markers(base_markers)

            if base_domain != domain:
                print(f'WARNING @ 0: {gid} has domain {base_domain} but metadata says {domain}')

            f.write(f'{gid}\t0\t{base_msa}\n')

            for cur_pct, contigs_to_omit in d_pct_to_contigs_to_remove.items():
                cur_results_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th, contigs_to_omit)
                cur_msa = get_msa_from_hits(cur_results_markers, marker_dir=marker_dir)
                cur_domain = get_genome_domain_from_markers(cur_results_markers)

                if cur_domain != domain:
                    print(f'WARNING @ {cur_pct}: {gid} has domain {cur_domain} but metadata says {domain}')

                f.write(f'{gid}\t{cur_pct}\t{cur_msa}\n')

        copy_file(path_tmp, path_out, checksum=True)

    return
