import multiprocessing as mp
import os
import tempfile
from collections import defaultdict
from typing import Tuple, List
import pandas as pd
from collections import defaultdict

from Bio import SeqIO
from magna.gunc import read_contig_assignments_tsv

from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.method.contig_removal import contigs_to_remove_from_gunc, get_taxonomy_by_majority_vote_gunc
from workflow.method.get_genome_domain_from_markers import get_genome_domain_from_markers
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid, get_all_marker_hits_for_gid
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.paths import get_gid_root
import os

from luigi import LocalTarget
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, PPLACER_R207_ARC_PATH, \
    DIR_OUT_PPLACER_ARC_NON_REPS, DIR_OUT_PPLACER_BAC_NON_REPS, PPLACER_R207_BAC_PATH, DIR_OUT_SENTINEL, \
    DIR_OUT_PPLACER_WGA, R207_BAC120_HMM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.marker.get_msa_for_failed_at_pct import GetMsaForFailedAtPct
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.model.pplacer import Pplacer
from workflow.tree.pplacer_non_reps_at_pct_arc_get_taxonomy import PplacerNonRepsAtPctArcGetTaxonomy
from workflow.tree.pplacer_non_reps_at_pct_bac_get_taxonomy import PplacerNonRepsAtPctBacGetTaxonomy
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_and_wait


class PplacerWgaOnChanged(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
            'arc_pplacer': PplacerNonRepsAtPctArcGetTaxonomy(),
            'bac_pplacer': PplacerNonRepsAtPctBacGetTaxonomy(),
        }

    def output(self):
        return LocalTarget(os.path.join(DIR_OUT_SENTINEL, self.fqn))

    def load_input(self):
        out = defaultdict(dict)
        df_bac = self.input()['arc_pplacer'].read() if not DEBUG else self.input()['arc_pplacer'].read_cached()
        df_arc = self.input()['bac_pplacer'].read() if not DEBUG else self.input()['bac_pplacer'].read_cached()
        df_merged = pd.concat([df_bac, df_arc], ignore_index=True)
        for row in df_merged.itertuples():
            out[row.gid][row.pct] = row.tax
        return out

    def run(self):
        log('Running whole genome alignment on those that changed genus', title=True)
        self.make_output_dirs()

        # DIR_OUT_PPLACER_WGA

        self.make_output_dirs()
        log('Loading metadata')
        df_meta = self.input()['meta'].read() if not DEBUG else self.input()['meta'].read_cached()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Merging dataframes')
        df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

        # log('Loading output from pplacer step')
        # d_gid_to_pct_tax = self.load_input()
        #
        # log('Loading gids that differ')
        # gids_that_differ = get_gids_that_differ(df_meta, d_gid_to_pct_tax)

        # return

        log('Manually inspect dataframe to select rows..')
        input_data = [
            {
                'gid': 'GCA_002170165.1',
                'to': [
                    'GCA_014381465.1',  # from family
                    'GCA_002724505.1',  # contaminant

                ],
                'algo': 'promer'
            },
            {
                'gid': 'GCA_002401495.1',
                'to': [
                    'GCA_017134395.1',  # from genus
                    'GCF_000009365.1',  # possible other contaminant
                    'GCF_002310815.1',  # possible other contaminant
                    'GCF_001613795.1',  # contaminant
                ],
                'algo': 'promer'
            },
            {
                'gid': 'GCA_014385135.2',
                'to': [
                    'GCF_014290175.1',  # species rep
                    'GCF_900291955.1',  # contaminant (hq)
                ]
            },
            {
                'gid': 'GCA_017515185.1',
                'to': [
                    'GCA_017538895.1',  # sp rep
                    'GCA_017509245.1',  # contaminant
                ],
                'algo': 'promer'
            },
            {
                'gid': 'GCA_900765805.1',
                'to': [
                    'GCF_013009555.1',  # sp rep
                    'GCF_000025985.1',  # contaminant
                ]

            },
            {
                'gid': 'GCA_905214645.1',
                'to': [
                    'GCF_014334015.1',  # sp rep
                    'GCF_000154105.1',  # contaminant
                ]
            }
        ]

        input_data_tmp = [
            {
                'gid': 'GCA_002170165.1',
                'to': ['GCA_002724505.1', 'GCA_014381465.1', 'GCA_012960835.1', 'GCA_002689995.1', 'GCA_002167775.1'],
            },
        ]

        input_data_tmp = [
            {
                'gid': 'GCA_000153745.1',
                'to': ['GCA_900197625.1', 'GCA_905182285.1', 'GCA_002292215.1', 'GCF_002807665.1'],
                'algo': 'promer'
            },
        ]


        log('Generating mummer2circos commands...')
        for job in input_data_tmp:
            cur_gid = job['gid']

            row = df_merged.loc[cur_gid]
            d_pct_to_contigs_to_remove = get_contig_contam_prop(cur_gid, row['source'],
                                                                row['domain'], row['taxonomic_level'])
            generate_command(f'{cur_gid}_bioproject',
                             cur_gid, job['to'],
                             d_pct_to_contigs_to_remove, job.get('algo'))

        return


def get_marker_locations(gid):
    out = list()

    # Read the FASTA file
    gid_root = get_gid_root(gid)

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_faa = SeqIO.to_dict(SeqIO.parse(path_faa, "fasta"))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    d_marker_to_pos, marker_hits = get_all_marker_hits_for_gid(d_faa, pfam_th, tigr_th)

    for marker in R207_BAC120_HMM:
        lst_pos = d_marker_to_pos.get(marker)
        if not lst_pos:
            continue

        if marker in marker_hits['unq'] or marker in marker_hits['muq']:
            label = f'{marker}'
        elif marker in marker_hits['mul']:
            label = f'{marker}'
        else:
            raise Exception('???')

        for i, (gene_id, from_pos, to_pos) in enumerate(lst_pos):
            contig_id = gene_id[:gene_id.rindex('_')]
            cur_label = label
            if marker in marker_hits['mul']:
                cur_label += f'-{i + 1}/{len(lst_pos)}'
            out.append((contig_id, from_pos, to_pos, cur_label))

    return out


def generate_command(dir_name, gid, to_gids, d_pct_to_contigs_to_remove, algo=None):
    out_dir = os.path.join(DIR_OUT_PPLACER_WGA, dir_name)
    os.makedirs(out_dir, exist_ok=True)

    d_contig_to_pcts = defaultdict(set)
    for pct, contigs in d_pct_to_contigs_to_remove.items():
        for contig in contigs:
            d_contig_to_pcts[contig].add(pct)

    marker_locations = get_marker_locations(gid)


    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')
    d_fna = read_fasta(gid_fna)

    # Re-order the contigs in the fasta file to match the order in which they're removed
    path_fna_reordered = os.path.join(out_dir, f'{gid}_reordered.fna')

    contig_ordering = list()
    for pct, contigs_to_remove in d_pct_to_contigs_to_remove.items():
        for contig in sorted(contigs_to_remove):
            if contig not in contig_ordering:
                contig_ordering.append(contig)
    for contig in d_fna:
        if contig not in contig_ordering:
            contig_ordering.append(contig)

    with open(path_fna_reordered, 'w') as f:
        for contig in contig_ordering:
            f.write(f'>{contig}\n{d_fna[contig]}\n')


    path_labels = os.path.join(out_dir, 'labels.tsv')
    with open(path_labels, 'w') as f:
        for row in marker_locations:
            f.write('\t'.join(map(str, row)) + '\n')

        # for contig, pcts in d_contig_to_pcts.items():
        #     contig_middle = len(d_fna[contig]) // 2
        #     pct_str = ','.join(map(str, sorted(pcts)))
        #     f.write(f'{contig}\t{contig_middle-1}\t{contig_middle+1}\t{pct_str}\n')

    cmd = [
        'mummer2circos',
        '-l',
        '-lf',
        path_labels,
        '-r',
        path_fna_reordered,
        '-q',
    ]

    download_path = os.path.join(out_dir, 'download.sh')
    to_download_cmds = list()
    for q_gid in tqdm(to_gids):
        q_gid_path = os.path.join(get_gid_root(q_gid), f'{q_gid}.fna')
        if not os.path.exists(q_gid_path):
            if q_gid.startswith('GCA'):
                q_gid_new = f'GCF{q_gid[3:]}'
            else:
                q_gid_new = f'GCA{q_gid[3:]}'
            q_gid_path = os.path.join(get_gid_root(q_gid_new), f'{q_gid_new}.fna')
        if not os.path.exists(q_gid_path):
            new_root = get_gid_root(q_gid, '/srv/home/uqamussi/projects/gunc-chimeras/output/genomes/aux')
            new_fna = os.path.join(new_root, f'{q_gid}.fna')
            to_download_cmds.append(f'magna ncbi download {q_gid} {new_fna}.gz && gunzip {new_fna}.gz')
            q_gid_path = new_fna
        cmd.append(q_gid_path)

    if len(to_download_cmds) > 0:
        with open(download_path, 'w') as f:
            f.write("\n".join(to_download_cmds) + '\n')

    if algo is not None:
        cmd.append('-a')
        cmd.append(algo)
    print(f'Working directory: {out_dir}')

    path_cmd = os.path.join(out_dir, 'cmd.sh')
    with open(path_cmd, 'w') as f:
        f.write(' '.join(cmd) + '\n')

    # print(' '.join(cmd) + '\n')

    return


def get_gids_that_differ(df_meta, d_gid_to_pct_tax):
    unq_gids = set()
    rows = list()
    for gid, d_pct_tax in sorted(d_gid_to_pct_tax.items()):
        expected_tax = df_meta.loc[gid]['gtdb_taxonomy'].split(';')
        for pct, tax in sorted(d_pct_tax.items()):
            cur_tax = tax.split(';')
            if expected_tax[0:6] != cur_tax[0:6]:
                print(f'{gid} @ {pct}')
                print(f'Expected: {";".join(expected_tax)}')
                print(f'Result  : {";".join(cur_tax)}')
                print()
                unq_gids.add(gid)
                rows.append({
                    'gid': gid,
                    'pct': pct,
                    'expected_tax': ';'.join(expected_tax),
                    'pct_tax': ';'.join(cur_tax)
                })

    print('Gids that were affected:')
    print('\n'.join(sorted(unq_gids)))
    df = pd.DataFrame(rows)
    return df.sort_values(['gid', 'pct'], ignore_index=True)


def get_contig_contam_prop(gid, gunc_source, domain, max_css):
    # Read the FASTA file
    gid_root = get_gid_root(gid)
    gid_fna = os.path.join(gid_root, f'{gid}.fna')

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
    return d_pct_to_contigs_to_remove
