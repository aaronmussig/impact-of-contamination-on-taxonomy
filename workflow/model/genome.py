import os
import random
import re
import tempfile
from collections import Counter, defaultdict
from functools import cached_property, lru_cache
from typing import Dict, Optional, List, FrozenSet, Union, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from gunc.get_scores import create_base_data as gunc_create_base_data
from gunc.get_scores import get_abundant_lineages_cutoff as gunc_get_abundant_lineages_cutoff
from gunc.get_scores import get_stats as gunc_get_stats
from gunc.get_scores import read_diamond_output as gunc_read_diamond_output
from magna.gunc import read_contig_assignments_tsv

from workflow.config import R207_MARKERS, R207_AR53_HMM, R207_BAC120_HMM, PATH_R207_AR53_MASK, PATH_R207_BAC120_MASK, \
    R207_MARKER_LENGTHS
from workflow.exception import GenomeNoMarkersFound
from workflow.method.contig_removal import get_contig_metadata
from workflow.method.get_msa_from_hits import align_marker
from workflow.model.gunc_model import GuncRefDb, GUNC_RANKS
from workflow.model.taxonomy import TaxDomain, Taxonomy
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile, Hit
from workflow.util.paths import get_gid_root


class Genome:
    """
    This class wraps all things to do with a single genome object.
    Computationally inefficient as things may need to be re-calculated multiple times,
    but much easier to keep track of what is going on.
    """

    # __slots__ = ('gid', 'root', 'memory')

    def __init__(self, gid: str):
        # JobLib requires that no objects in self change once the class
        # has been instantiated, or it renders the cache invalid!
        self.gid = gid
        self.root: str = get_gid_root(gid)

        # # Create the Joblib backend
        # self.memory = Memory(location=get_gid_root(gid, DIR_CACHE_GENOME), verbose=0)
        #
        # # Set the functions to be cached
        # self.get_merged_top_hit = self.memory.cache(self.get_merged_top_hit, ignore=['self'])
        # self.get_marker_hits = self.memory.cache(self.get_marker_hits, ignore=['self'])
        # self.get_domain_from_markers = self.memory.cache(self.get_domain_from_markers, ignore=['self'])
        # self.get_aligned_unq_domain_markers = self.memory.cache(self.get_aligned_unq_domain_markers, ignore=['self'])
        # self.get_domain_msa = self.memory.cache(self.get_domain_msa, ignore=['self'])
        # self.get_gunc_max_css_inferred_taxonomy = self.memory.cache(self.get_gunc_max_css_inferred_taxonomy, ignore=['self'])
        # self.get_gunc_max_css_contig_ranking = self.memory.cache(self.get_gunc_max_css_contig_ranking, ignore=['self'])
        # self.get_gunc_contigs_where_removed_equals_x_pct_genome_removed = self.memory.cache(self.get_gunc_contigs_where_removed_equals_x_pct_genome_removed, ignore=['self'])
        # self.get_marker_congruence = self.memory.cache(self.get_marker_congruence, ignore=['self'])

    def __repr__(self):
        return self.gid

    @cached_property
    def d_fna(self) -> Dict[str, str]:
        """Return the contig to DNA sequences for this genome."""
        out = dict()
        path_fna: str = os.path.join(self.root, f'{self.gid}.fna')
        with open(path_fna) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                out[record.id] = str(record.seq)
        return out

    @cached_property
    def d_faa(self) -> Dict[str, str]:
        """Return the contig to gene sequences for this genome."""
        out = dict()
        path_faa = os.path.join(self.root, 'prodigal', f'{self.gid}.faa')
        with open(path_faa) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.seq.endswith('*'):
                    out[record.id] = str(record.seq[:-1])
                else:
                    out[record.id] = str(record.seq)
        return out

    @cached_property
    def pfam_top_hit(self) -> TopHitPfamFile:
        """Return the PFAM TopHit file."""
        path = os.path.join(self.root, 'prodigal', f'{self.gid}_pfam_tophit.tsv')
        file = TopHitPfamFile(path)
        file.read()
        return file

    @cached_property
    def tigrfam_top_hit(self) -> TopHitTigrFile:
        """Return the TIGRFAM TopHit file."""
        path = os.path.join(self.root, 'prodigal', f'{self.gid}_tigrfam_tophit.tsv')
        file = TopHitTigrFile(path)
        file.read()
        return file

    @lru_cache
    def get_merged_top_hit(self, omit_contigs: Optional[FrozenSet[str]] = None) -> Dict[str, List[Hit]]:
        """Merge the top hit files, omitting contigs if provided."""
        if omit_contigs is None:
            omit_contigs = set()
        out = dict()
        for cur_tophit in (self.pfam_top_hit, self.tigrfam_top_hit):
            for gene_id, cur_hit in cur_tophit.iter_hits():
                contig_id = gene_id[:gene_id.rindex('_')]
                if contig_id not in omit_contigs:
                    if cur_hit.hmm_id not in out:
                        out[cur_hit.hmm_id] = list()
                    out[cur_hit.hmm_id].append(cur_hit)
        return out

    @lru_cache
    def get_marker_hits(self, omit_contigs: Optional[FrozenSet[str]] = None):
        """Calculate the unique, multiple, missing, duplicate marker hits. Returns the R207 markers only."""

        # Pointers to unique, multiple hit, multiple-unique, missing markers.
        cur_unq = dict()
        cur_mul = dict()
        cur_muq = dict()
        cur_mis = set()

        # Create a dictionary of marker names -> Hits
        d_hmm_hits = self.get_merged_top_hit(omit_contigs=omit_contigs)

        # Foreach expected marker determine which category it falls into.
        for marker_id in R207_MARKERS:

            # Marker is missing.
            if marker_id not in d_hmm_hits:
                cur_mis.add(marker_id)

            # Multiple hits to the same marker.
            elif len(d_hmm_hits[marker_id]) > 1:

                # If sequences are the same, take the most significant hit
                unq_seqs = {self.d_faa[x.gene_id] for x in d_hmm_hits[marker_id]}
                if len(unq_seqs) == 1:
                    cur_top_hit = sorted(d_hmm_hits[marker_id], reverse=True)[0]
                    cur_muq[marker_id] = {'hit': cur_top_hit, 'seq': self.d_faa[cur_top_hit.gene_id]}

                # Marker maps to multiple genes.
                else:
                    if marker_id not in cur_mul:
                        cur_mul[marker_id] = list()
                    for x_i, x in enumerate(d_hmm_hits[marker_id]):
                        cur_mul[marker_id].append({'hit': x, 'seq': self.d_faa[x.gene_id]})

            # This was a unique hit.
            else:
                cur_hit = d_hmm_hits[marker_id][0]
                cur_unq[marker_id] = {'hit': cur_hit, 'seq': self.d_faa[cur_hit.gene_id]}

        # invert the dict
        out_markers = {
            'unq': cur_unq,
            'mul': cur_mul,
            'muq': cur_muq,
            'mis': cur_mis,
        }
        return out_markers

    @lru_cache
    def get_domain_from_markers(self, omit_contigs: Optional[FrozenSet[str]] = None) -> Dict[str, Union[str, float]]:
        """Return the domain of the genome based on the markers."""
        markers = self.get_marker_hits(omit_contigs=omit_contigs)
        single_copy_hits = {**markers['muq'], **markers['unq']}

        arc_count = len(set(single_copy_hits).intersection(set(R207_AR53_HMM)))
        bac_count = len(set(single_copy_hits).intersection(set(R207_BAC120_HMM)))

        arc_aa_pct = arc_count / len(R207_AR53_HMM) * 100
        bac_aa_pct = bac_count / len(R207_BAC120_HMM) * 100

        if arc_count == 0 and bac_count == 0:
            raise GenomeNoMarkersFound(f'No markers found for: {self.gid}')

        domain = TaxDomain.BACTERIA if bac_aa_pct >= arc_aa_pct else TaxDomain.ARCHAEA
        return {
            'domain': domain,
            'arc_aa_pct': arc_aa_pct,
            'bac_aa_pct': bac_aa_pct,
        }

    @lru_cache
    def get_aligned_unq_domain_markers(self, omit_contigs: Optional[FrozenSet[str]] = None) -> Dict[str, str]:
        """Return the aligned markers for this genome."""

        # Get the domain
        domain = self.get_domain_from_markers(omit_contigs=omit_contigs)['domain']
        if domain is TaxDomain.ARCHAEA:
            marker_ids = R207_AR53_HMM
        elif domain is TaxDomain.BACTERIA:
            marker_ids = R207_BAC120_HMM
        else:
            raise ValueError(f'Unknown domain: {domain}')

        # Get the unique marker hits.
        marker_hits = self.get_marker_hits(omit_contigs=omit_contigs)
        single_copy_hits = {**marker_hits['muq'], **marker_hits['unq']}
        markers_to_align = set(single_copy_hits).intersection(set(marker_ids))

        # Align the markers
        d_marker_to_aln = dict()
        with tempfile.TemporaryDirectory(prefix='gunc-chim-msa') as tmp_dir:
            for marker in sorted(markers_to_align):
                hit = single_copy_hits[marker]

                path = os.path.join(tmp_dir, f'{marker}.faa')
                with open(path, 'w') as f:
                    f.write(f'>{marker}\n{hit["seq"]}\n')

                aligned_markers = align_marker(marker, path)
                d_marker_to_aln[marker] = aligned_markers[marker]

        # Fill in the missing markers with gaps
        out = dict()
        for marker in marker_ids:
            if marker not in d_marker_to_aln:
                out[marker] = '-' * R207_MARKER_LENGTHS[marker]
            else:
                out[marker] = d_marker_to_aln[marker]

        return out

    def get_aligned_unq_domain_markers_from_markers(self, d_markers: dict, masked=False) -> Tuple[Dict[str, str], str]:
        """Return the aligned markers for this genome from a specific input."""

        # Get the domain
        domain = self.get_domain_from_markers()['domain']
        if domain is TaxDomain.ARCHAEA:
            marker_ids = R207_AR53_HMM
        elif domain is TaxDomain.BACTERIA:
            marker_ids = R207_BAC120_HMM
        else:
            raise ValueError(f'Unknown domain: {domain}')

        # Align the markers
        d_marker_to_aln = dict()
        with tempfile.TemporaryDirectory(prefix='gunc-chim-msa') as tmp_dir:
            for marker, hit in sorted(d_markers.items()):
                path = os.path.join(tmp_dir, f'{marker}.faa')
                with open(path, 'w') as f:
                    f.write(f'>{marker}\n{hit["seq"]}\n')

                aligned_markers = align_marker(marker, path)
                d_marker_to_aln[marker] = aligned_markers[marker]

        # Fill in the missing markers with gaps
        out = dict()
        for marker in marker_ids:
            if marker not in d_marker_to_aln:
                out[marker] = '-' * R207_MARKER_LENGTHS[marker]
            else:
                out[marker] = d_marker_to_aln[marker]

        # Load the mask (if masking)
        if masked:
            if domain is TaxDomain.ARCHAEA:
                path_mask = PATH_R207_AR53_MASK
            elif domain is TaxDomain.BACTERIA:
                path_mask = PATH_R207_BAC120_MASK
            else:
                raise ValueError(f'Unknown domain: {domain}')
            with open(path_mask) as f:
                mask = [x == '1' for x in f.read().strip()]
        else:
            mask = None

        # Create the MSA
        msa = ''.join([v for k, v in sorted(out.items(), key=lambda x: x[0])])

        # Trim if required
        if mask is not None:
            msa = ''.join([x for x, y in zip(msa, mask) if y])

        return out, msa

    @lru_cache
    def get_domain_msa(self, masked=True, omit_contigs: Optional[FrozenSet[str]] = None) -> str:
        """Return the multiple sequence alignment for the inferred domain (by marker count)."""

        # Load the mask (if masking)
        if masked:
            domain = self.get_domain_from_markers(omit_contigs=omit_contigs)['domain']
            if domain is TaxDomain.ARCHAEA:
                path_mask = PATH_R207_AR53_MASK
            elif domain is TaxDomain.BACTERIA:
                path_mask = PATH_R207_BAC120_MASK
            else:
                raise ValueError(f'Unknown domain: {domain}')
            with open(path_mask) as f:
                mask = [x == '1' for x in f.read().strip()]
        else:
            mask = None

        # Get the unique marker hits.
        d_marker_to_aln = self.get_aligned_unq_domain_markers(omit_contigs=omit_contigs)

        # Create the MSA
        msa = ''.join([v for k, v in sorted(d_marker_to_aln.items(), key=lambda x: x[0])])

        # Trim if required
        if mask is not None:
            msa = ''.join([x for x, y in zip(msa, mask) if y])
        return msa

    def get_random_contigs_where_removed_equals_x_pct_genome_removed(self, pct: float) -> Tuple[List[str], float]:
        """Return the list of contigs, that when removed equals x pct of the genome removed."""

        total_n_bases = sum([len(x) for x in self.d_fna.values()])

        # Shuffle the contigs
        contigs = list(self.d_fna.keys())
        random.shuffle(contigs)

        # Remove the contigs in order until we hit the pct
        out = list()
        nt_removed = 0
        for contig in contigs:
            contig_len = len(self.d_fna[contig])

            new_pct_removed = 100 * (nt_removed + contig_len) / total_n_bases

            if new_pct_removed >= 100:
                break

            nt_removed += contig_len
            out.append(contig)

            if new_pct_removed >= pct:
                break

        pct_removed = 100 * nt_removed / total_n_bases
        return out, pct_removed

    @cached_property
    def df_gunc_contig_assignment_r95(self) -> pd.DataFrame:
        """Return the gunc contig assignment table."""
        path = os.path.join(self.root, 'gunc_r95', 'gunc_output', f'{self.gid}.contig_assignments.tsv')
        return read_contig_assignments_tsv(path)

    @cached_property
    def df_gunc_contig_assignment_pro(self) -> pd.DataFrame:
        """Return the gunc contig assignment table."""
        path = os.path.join(self.root, 'gunc_pro', 'gunc_output', f'{self.gid}.contig_assignments.tsv')
        return read_contig_assignments_tsv(path)

    @cached_property
    def path_gunc_diamond_r95(self) -> str:
        """Return the path to the gunc diamond output."""
        return os.path.join(self.root, 'gunc_r95', 'diamond_output', f'{self.gid}.diamond.gtdb_95.out')

    @cached_property
    def path_gunc_diamond_pro(self) -> str:
        """Return the path to the gunc diamond output."""
        return os.path.join(self.root, 'gunc_pro', 'diamond_output', f'{self.gid}.diamond.progenomes_2.1.out')

    @lru_cache
    def get_gunc_max_css_inferred_taxonomy(self, max_css: GUNC_RANKS, source: GuncRefDb) -> Taxonomy:
        """Return the gunc inferred taxonomy (only up to the MaxCSS). Replicates GUNC logic."""

        # Load the appropriate diamond file
        if source is GuncRefDb.GTDB:
            diamond_file_path = self.path_gunc_diamond_r95
            gunc_db_str = 'gtdb_95'
        elif source is GuncRefDb.PRO:
            diamond_file_path = self.path_gunc_diamond_pro
            gunc_db_str = 'progenomes_2.1'
        else:
            raise ValueError(f'Unknown source: {source}')

        # This repeats the method used in the GUNC program
        diamond_df = gunc_read_diamond_output(diamond_file_path)
        base_data = gunc_create_base_data(diamond_df, gunc_db_str)

        genes_mapped, contig_count = gunc_get_stats(diamond_df)
        abundant_lineages_cutoff = gunc_get_abundant_lineages_cutoff(False, genes_mapped)

        tax_data = base_data.groupby(max_css).filter(
            lambda x: len(x) > abundant_lineages_cutoff
        )

        df_taxstrings = tax_data[list(GUNC_RANKS[:GUNC_RANKS.index(max_css) + 1])]
        taxstrings = [';'.join(x) for x in df_taxstrings.values]
        taxstrings_count = Counter(taxstrings)

        # Can't reach a consensus if there are multiple most common
        if len(taxstrings_count) > 1 and len(set(x[1] for x in taxstrings_count.most_common(2))) == 1:
            raise Exception(f'{self.gid} has multiple most common taxstrings')

        # Get the most common taxstring
        taxonomy = taxstrings_count.most_common(1)[0][0]
        return Taxonomy(*taxonomy.split(';'))

    def get_gunc_max_css_contig_assignments(self, max_css: GUNC_RANKS, source: GuncRefDb):
        # Get the taxonomy, and the specific taxon at this max css level
        taxonomy = self.get_gunc_max_css_inferred_taxonomy(max_css=max_css, source=source)
        taxon = taxonomy.from_gunc_rank(max_css)

        # Load the appropriate contig assignment file
        if source is GuncRefDb.GTDB:
            df_contig_assign = self.df_gunc_contig_assignment_r95
        elif source is GuncRefDb.PRO:
            df_contig_assign = self.df_gunc_contig_assignment_pro
        else:
            raise ValueError(f'Unknown source: {source}')

        d_contig_to_assignments = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
        for row in df_contig_assign.itertuples():
            d_contig_to_assignments[row.contig][row.tax_level][row.assignment] += row.count_of_genes_assigned

        out = dict()
        for contig, d_tax_level_to_assignment in d_contig_to_assignments.items():
            if len(d_tax_level_to_assignment['genus']) > 0:
                out[contig] = sorted(d_tax_level_to_assignment['genus'].items(), key=lambda x: -x[1])[0][0]
        return out

    @lru_cache
    def get_gunc_max_css_marker_ranking(self, max_css: GUNC_RANKS, source: GuncRefDb):
        """Return the contig ranking for the given MaxCSS."""

        # Get the taxonomy, and the specific taxon at this max css level
        taxonomy = self.get_gunc_max_css_inferred_taxonomy(max_css=max_css, source=source)

        # Load the appropriate diamond file
        if source is GuncRefDb.GTDB:
            diamond_file_path = self.path_gunc_diamond_r95
            gunc_db_str = 'gtdb_95'
        elif source is GuncRefDb.PRO:
            diamond_file_path = self.path_gunc_diamond_pro
            gunc_db_str = 'progenomes_2.1'
        else:
            raise ValueError(f'Unknown source: {source}')

        # This repeats the method used in the GUNC program
        diamond_df = gunc_read_diamond_output(diamond_file_path)
        base_data = gunc_create_base_data(diamond_df, gunc_db_str)

        genes_mapped, contig_count = gunc_get_stats(diamond_df)
        abundant_lineages_cutoff = gunc_get_abundant_lineages_cutoff(False, genes_mapped)

        tax_data = base_data.groupby(max_css).filter(
            lambda x: len(x) > abundant_lineages_cutoff
        )

        # Load the marker hits and flatten into a list
        marker_hits = self.get_marker_hits()
        marker_hits_list = list()
        for _, d_hit in {**marker_hits['unq'], **marker_hits['muq']}.items():
            marker_hits_list.append(d_hit)
        for _, lst_d_hit in marker_hits['mul'].items():
            marker_hits_list.extend(lst_d_hit)

        # Load the contig ranking for additional ranking information
        contig_ranking = self.get_gunc_max_css_contig_ranking(max_css, source)

        # Iterate over marker hits to see if they're correct
        rows = list()
        for i, d_hit in enumerate(marker_hits_list):
            cur_gene_id = d_hit['hit'].gene_id
            cur_contig = cur_gene_id[:cur_gene_id.rindex('_')]

            cur_contig_ranking = contig_ranking[contig_ranking['contig'] == cur_contig]
            if len(cur_contig_ranking) != 1:
                raise Exception('?')
            cur_contig_ranking = cur_contig_ranking.iloc[0]

            cur_diamond_hit = tax_data[tax_data['query'] == cur_gene_id]
            if len(cur_diamond_hit) > 1:
                raise Exception('?')
            elif len(cur_diamond_hit) == 0:
                cur_diamond_hit = {
                    'kingdom': 'NA',
                    'phylum': 'NA',
                    'class': 'NA',
                    'order': 'NA',
                    'family': 'NA',
                    'genus': 'NA',
                    'id': '-1',
                }
            else:
                cur_diamond_hit = cur_diamond_hit.iloc[0]

            cur_row = {
                'idx': i,
                'marker': d_hit['hit'].hmm_id,
                'e_val': d_hit['hit'].e_val,
                'bit_score': d_hit['hit'].bit_score,
                'gene_id': cur_gene_id,
                'contig': cur_contig,

                'g_kingdom': 1 if taxonomy.d == cur_diamond_hit['kingdom'] else 0,
                'g_phylum': 1 if taxonomy.p == cur_diamond_hit['phylum'] else 0,
                'g_class': 1 if taxonomy.c == cur_diamond_hit['class'] else 0,
                'g_order': 1 if taxonomy.o == cur_diamond_hit['order'] else 0,
                'g_family': 1 if taxonomy.f == cur_diamond_hit['family'] else 0,
                'g_genus': 1 if taxonomy.g == cur_diamond_hit['genus'] else 0,
                'identity_diamond': cur_diamond_hit['id'],

            }

            d_other_to_add = cur_contig_ranking.to_dict()
            del d_other_to_add['order']
            cur_row = {**cur_row, **d_other_to_add}

            # Save the result
            rows.append(cur_row)

        # Create temporary output dataframe
        df_out = pd.DataFrame(rows)

        # Get the higher ranks to the max css level
        gunc_ranks_to_max_css = GUNC_RANKS[0:GUNC_RANKS.index(max_css) + 1]

        # Order by the % of congruence for each contig
        columns_to_sort = list()
        columns_to_sort.extend([f'pct_{x}' for x in reversed(gunc_ranks_to_max_css)])
        columns_to_sort.extend([f'g_{x}' for x in reversed(gunc_ranks_to_max_css)])
        columns_to_sort.extend(['length', 'identity_diamond', 'bit_score', 'contig', 'marker'])

        # Order the results
        df_out.sort_values(
            by=columns_to_sort,
            ascending=True,
            inplace=True,
            ignore_index=True
        )

        return df_out, marker_hits_list

    @lru_cache
    def get_gunc_max_css_contig_ranking(self, max_css: GUNC_RANKS, source: GuncRefDb):
        """Return the contig ranking for the given MaxCSS."""

        # Get the taxonomy, and the specific taxon at this max css level
        taxonomy = self.get_gunc_max_css_inferred_taxonomy(max_css=max_css, source=source)
        taxon = taxonomy.from_gunc_rank(max_css)

        # Load the appropriate contig assignment file
        if source is GuncRefDb.GTDB:
            df_contig_assign = self.df_gunc_contig_assignment_r95
        elif source is GuncRefDb.PRO:
            df_contig_assign = self.df_gunc_contig_assignment_pro
        else:
            raise ValueError(f'Unknown source: {source}')

        # Get the higher ranks to the max css level
        gunc_ranks_to_max_css = GUNC_RANKS[0:GUNC_RANKS.index(max_css) + 1]

        # Iterate over each higher rank and determine the % mapping to that taxon
        # this has the benefit of tie-breaking 0% maps better, for more similar taxa
        d_contig_to_taxon_counts = defaultdict(lambda: defaultdict(lambda: {'correct': 0, 'total': 0}))
        for gunc_rank in gunc_ranks_to_max_css:
            expected_taxon = taxonomy.from_gunc_rank(gunc_rank)

            # Subset the dataframe to only this rank
            df_contig_assign_subset = df_contig_assign[df_contig_assign['tax_level'] == gunc_rank]

            # Count the number of correct vs incorrect hits
            for row in df_contig_assign_subset.itertuples():
                cnt_assigned = int(row.count_of_genes_assigned)
                if row.assignment == expected_taxon:
                    d_contig_to_taxon_counts[row.contig][gunc_rank]['correct'] += cnt_assigned
                d_contig_to_taxon_counts[row.contig][gunc_rank]['total'] += cnt_assigned

        # Create the ranking dataframe (still need to get some contig metadata)
        out = list()

        for contig, contig_metadata in get_contig_metadata(self.d_fna).items():
            d_taxon_counts = d_contig_to_taxon_counts[contig]

            cur_row = {
                'gid': self.gid,
                'contig': contig,
                'length': contig_metadata['length'],
                'gc': contig_metadata['gc'],
            }

            if d_taxon_counts[max_css]['total'] == 0:
                cur_row['pct_correct'] = 0
            else:
                cur_row['pct_correct'] = 100 * d_taxon_counts[max_css]['correct'] / d_taxon_counts[max_css]['total']

            # Iterate over each higher rank to determine the mapping %
            for gunc_rank in reversed(gunc_ranks_to_max_css):
                d_counts = d_taxon_counts[gunc_rank]
                if d_counts['total'] == 0:
                    cur_row[f'pct_{gunc_rank}'] = 0
                else:
                    cur_row[f'pct_{gunc_rank}'] = 100 * d_counts['correct'] / d_counts['total']
            out.append(cur_row)

        # Create the dataframe
        df_out = pd.DataFrame(out)

        # Order by the % of congruence for each contig
        columns_to_sort = [f'pct_{x}' for x in reversed(gunc_ranks_to_max_css)]
        columns_to_sort.append('length')
        columns_to_sort.append('contig')
        df_out.sort_values(by=columns_to_sort, inplace=True, ignore_index=True)

        # Create the ordering
        df_out['order'] = np.arange(0, len(df_out))
        return df_out

    @lru_cache
    def get_gunc_contigs_where_removed_equals_x_pct_genome_removed(
            self,
            pct: float,
            max_congruence: float,
            max_css: GUNC_RANKS,
            source: GuncRefDb
    ):

        # Get the GUNC contig ranking
        df_contig_ranking = self.get_gunc_max_css_contig_ranking(max_css=max_css, source=source)

        # Get the total number of bases in the genome
        total_n_bases = sum([len(x) for x in self.d_fna.values()])

        # Remove the contigs in order until we hit the pct
        out = list()
        nt_removed = 0
        for row in df_contig_ranking.itertuples():
            contig = str(row.contig)
            contig_len = int(row.length)

            # If we're removing contigs that map well, stop
            if row.pct_correct >= max_congruence:
                break

            # Calculate what percentage we will be at by removing this contig
            new_pct_removed = 100 * (nt_removed + contig_len) / total_n_bases

            # Don't remove all contigs
            if new_pct_removed >= 100:
                break

            # Add the contig so we can at least get to the target
            nt_removed += contig_len
            out.append(contig)

            # Check if we have exceeded the goal
            if new_pct_removed >= pct:
                break

        pct_removed = 100 * nt_removed / total_n_bases
        return out, pct_removed

    @lru_cache
    def get_marker_congruence(self, max_css: GUNC_RANKS, source: GuncRefDb,
                              omit_contigs: Optional[FrozenSet[str]] = None):
        """Return the marker congruence for the given MaxCSS."""

        marker_hits = self.get_marker_hits(omit_contigs)
        contig_ranking = self.get_gunc_max_css_contig_ranking(max_css, source)

        unq_hits = {**marker_hits['unq'], **marker_hits['muq']}

        d_contig_to_pct_correct = contig_ranking[['contig', 'pct_correct']].set_index('contig').to_dict()['pct_correct']

        out = dict()
        for marker, hit in sorted(unq_hits.items()):
            gene_id = hit['hit'].gene_id
            contig_id = gene_id[:gene_id.rindex('_')]
            out[marker] = d_contig_to_pct_correct[contig_id]

        return out

    @lru_cache
    def get_unq_markers_present_at_pct_removed(self, max_css: GUNC_RANKS, source: GuncRefDb, pct: float):
        """Return all markers ranked by their congruence value / position on the contig.
        Note: pct is the percent of markers to keep
        """

        # Get the taxonomy, and the specific taxon at this max css level
        taxonomy = self.get_gunc_max_css_inferred_taxonomy(max_css=max_css, source=source)

        # Get the marker set to use
        domain = self.get_domain_from_markers()['domain']
        if domain is TaxDomain.ARCHAEA:
            marker_ids = R207_AR53_HMM
        elif domain is TaxDomain.BACTERIA:
            marker_ids = R207_BAC120_HMM
        else:
            raise ValueError(f'Unknown domain: {domain}')
        marker_ids = frozenset(marker_ids)

        # Quick sanity check, the domain should match the majority vote
        if source is GuncRefDb.PRO:
            if taxonomy.d == '2157 Archaea':
                domain_vote = TaxDomain.ARCHAEA
            else:
                domain_vote = TaxDomain.BACTERIA
        else:
            if taxonomy.d == 'd__Archaea':
                domain_vote = TaxDomain.ARCHAEA
            else:
                domain_vote = TaxDomain.BACTERIA
        if domain_vote is not domain:
            raise Exception('??')

        # Calculate how many markers we start with (i.e. unique / multi-unique)
        markers_expected = set({**self.get_marker_hits()['unq'], **self.get_marker_hits()['muq']}.keys()).intersection(
            marker_ids)
        n_markers_expected = len(markers_expected)
        n_markers_to_keep = int(np.floor(n_markers_expected * (pct / 100)))

        # Get the marker congruency ranking
        df_ranking, marker_hits = self.get_gunc_max_css_marker_ranking(max_css=max_css, source=source)

        # Iterate over the markers in reverse order to build a list of the most
        # congruent unique markers that are from the expected domain
        d_marker_to_idx_to_keep_a = dict()
        d_marker_to_idx_to_keep_b = dict()
        markers_to_ignore_a = set()
        markers_to_ignore_b = set()
        cur_marker_to_idx_to_keep = d_marker_to_idx_to_keep_a
        cur_markers_to_ignore = markers_to_ignore_a

        # Collect markers
        for row in reversed(list(df_ranking.itertuples())):

            # Check to see if we have the target number of markers present
            if len(d_marker_to_idx_to_keep_a) >= n_markers_to_keep:
                # Swap the variables we're saving to (the second half of the chimera)
                cur_marker_to_idx_to_keep = d_marker_to_idx_to_keep_b
                cur_markers_to_ignore = markers_to_ignore_b

            # These are markers that are not within the domain markers
            if row.marker not in marker_ids:
                continue

            # These are markers that were marked as duplicate, and will be ignored
            if row.marker in cur_markers_to_ignore:
                continue

            # Check if it's duplicate/multi-unique
            if row.marker in cur_marker_to_idx_to_keep:

                # This marker already exists, check to see if it's unique or a duplicate
                idx_existing = cur_marker_to_idx_to_keep[row.marker]
                seq_existing = marker_hits[idx_existing]['seq']
                seq_current = marker_hits[row.idx]['seq']

                # This is a multi-unique hit, ignore it
                if seq_existing == seq_current:
                    continue

                # This is a duplicate hit, delete the old one
                else:
                    cur_markers_to_ignore.add(row.marker)
                    del cur_marker_to_idx_to_keep[row.marker]

            # No data on this hit yet, keep it
            else:
                cur_marker_to_idx_to_keep[row.marker] = row.idx

        # Make sure we got the same number back
        if len(d_marker_to_idx_to_keep_a) != n_markers_to_keep:
            raise Exception('?')

        # Output the markers
        d_out_a = dict()
        for marker, idx in d_marker_to_idx_to_keep_a.items():
            d_out_a[marker] = marker_hits[idx]
        d_out_b = dict()
        for marker, idx in d_marker_to_idx_to_keep_b.items():
            d_out_b[marker] = marker_hits[idx]
        return d_out_a, markers_expected, domain_vote, d_out_b

    def get_unq_random_markers_present_at_pct_removed(self, pct: float):
        """Return all markers randomly ranked / position on the contig."""

        # Get the marker set to use
        domain = self.get_domain_from_markers()['domain']
        if domain is TaxDomain.ARCHAEA:
            marker_ids = R207_AR53_HMM
        elif domain is TaxDomain.BACTERIA:
            marker_ids = R207_BAC120_HMM
        else:
            raise ValueError(f'Unknown domain: {domain}')
        marker_ids = frozenset(marker_ids)

        # Calculate how many markers we start with (i.e. unique / multi-unique)
        markers_expected = set({**self.get_marker_hits()['unq'], **self.get_marker_hits()['muq']}.keys()).intersection(
            marker_ids)
        n_markers_expected = len(markers_expected)
        n_markers_to_keep = int(np.floor(n_markers_expected * (pct / 100)))

        # Load the sequences for the markers
        marker_hits = self.get_marker_hits()
        d_gene_id_to_seq = dict()
        for marker_id, d_marker_hit in marker_hits['unq'].items():
            d_gene_id_to_seq[d_marker_hit['hit'].gene_id] = (d_marker_hit['seq'], marker_id)
        for marker_id, d_marker_hit in marker_hits['muq'].items():
            d_gene_id_to_seq[d_marker_hit['hit'].gene_id] = (d_marker_hit['seq'], marker_id)
        for marker_id, lst_marker_hit in marker_hits['mul'].items():
            for d_marker_hit in lst_marker_hit:
                d_gene_id_to_seq[d_marker_hit['hit'].gene_id] = (d_marker_hit['seq'], marker_id)

        # Randomly establish a contig ranking
        contigs = sorted(self.d_fna.keys())
        random.shuffle(contigs)

        # Load the marker locations
        d_contig_to_marker = defaultdict(list)
        for gene_id, (seq, marker_id) in d_gene_id_to_seq.items():
            cur_contig = gene_id[:gene_id.rindex('_')]
            d_contig_to_marker[cur_contig].append((marker_id, gene_id, seq))

        # Create the ordering
        lst_order = list()
        for contig in contigs:
            lst_markers = d_contig_to_marker[contig]
            if len(lst_markers) == 0:
                continue
            lst_markers = sorted(lst_markers, key=lambda x: int(x[1][x[1].rindex('_') + 1:]))
            lst_order.extend(lst_markers)

        # Iterate over the markers in reverse order to build a list of the most
        # congruent unique markers that are from the expected domain
        d_marker_to_idx_to_keep_a = dict()
        d_marker_to_idx_to_keep_b = dict()
        markers_to_ignore_a = set()
        markers_to_ignore_b = set()
        cur_marker_to_idx_to_keep = d_marker_to_idx_to_keep_a
        cur_markers_to_ignore = markers_to_ignore_a

        # Collect markers
        for idx, (marker_id, gene_id, seq) in enumerate(lst_order):

            # Check to see if we have the target number of markers present
            if len(d_marker_to_idx_to_keep_a) >= n_markers_to_keep:
                # Swap the variables we're saving to (the second half of the chimera)
                cur_marker_to_idx_to_keep = d_marker_to_idx_to_keep_b
                cur_markers_to_ignore = markers_to_ignore_b

            # These are markers that are not within the domain markers
            if marker_id not in marker_ids:
                continue

            # These are markers that were marked as duplicate, and will be ignored
            if marker_id in cur_markers_to_ignore:
                continue

            # Check if it's duplicate/multi-unique
            if marker_id in cur_marker_to_idx_to_keep:

                # This marker already exists, check to see if it's unique or a duplicate
                idx_existing = cur_marker_to_idx_to_keep[marker_id]
                seq_existing = lst_order[idx_existing][2]
                seq_current = lst_order[idx][2]

                # This is a multi-unique hit, ignore it
                if seq_existing == seq_current:
                    continue

                # This is a duplicate hit, delete the old one
                else:
                    cur_markers_to_ignore.add(marker_id)
                    del cur_marker_to_idx_to_keep[marker_id]

            # No data on this hit yet, keep it
            else:
                cur_marker_to_idx_to_keep[marker_id] = idx

        # Make sure we got the same number back
        if len(d_marker_to_idx_to_keep_a) != n_markers_to_keep:
            raise Exception('?')

        # Output the markers
        d_out_a = dict()
        for marker, idx in d_marker_to_idx_to_keep_a.items():
            d_out_a[marker] = {'seq': lst_order[idx][2]}
        d_out_b = dict()
        for marker, idx in d_marker_to_idx_to_keep_b.items():
            d_out_b[marker] = {'seq': lst_order[idx][2]}
        return d_out_a, markers_expected, domain, d_out_b

    @lru_cache
    def get_marker_loci(self) -> Dict[str, List[Tuple[str, int, int]]]:

        # Load the description of each gene call
        d_gene_to_description = dict()
        path_faa = os.path.join(self.root, 'prodigal', f'{self.gid}.faa')
        with open(path_faa) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                d_gene_to_description[record.id] = record.description

        re_pos = re.compile(r'.+? # (\d+) # (\d+) #.+')

        d_marker_to_position = defaultdict(list)

        # Create a dictionary of marker names -> Hits

        d_hmm_hits = self.get_merged_top_hit()

        # Foreach expected marker determine which category it falls into.
        for marker_id in R207_MARKERS:

            # Marker is missing.
            if marker_id not in d_hmm_hits:
                continue

            for hit in d_hmm_hits.get(marker_id, list()):
                description = d_gene_to_description[hit.gene_id]

                re_hit = re_pos.search(description)
                if re_hit is None:
                    raise Exception('??')
                pos_from, pos_to = re_hit.groups()
                pos_from, pos_to = int(pos_from), int(pos_to)
                d_marker_to_position[marker_id].append((hit.gene_id, pos_from, pos_to))

        return d_marker_to_position

    def get_gunc_max_css_marker_ranking_and_coords(self, max_css: GUNC_RANKS, source: GuncRefDb):

        # Get the marker ranking
        df_marker_ranking, _ = self.get_gunc_max_css_marker_ranking(max_css, source)

        # Obtain the marker loci (i.e. the start and end positions of each marker)
        d_marker_to_lst_loci = self.get_marker_loci()

        # Reverse the dictionary to obtain the gene id to coordinates
        d_gene_id_to_coord = dict()
        for lst_genes in d_marker_to_lst_loci.values():
            for gene_id, pos_from, pos_to in lst_genes:
                d_gene_id_to_coord[gene_id] = (pos_from, pos_to)

        # Append these to the marker ranking dataframe
        lst_pos_from, lst_pos_to = list(), list()
        for row in df_marker_ranking.itertuples():
            cur_pos_from, cur_pos_to = d_gene_id_to_coord[row.gene_id]
            lst_pos_from.append(cur_pos_from)
            lst_pos_to.append(cur_pos_to)
        df_marker_ranking['pos_from'] = lst_pos_from
        df_marker_ranking['pos_to'] = lst_pos_to
        df_marker_ranking['rank'] = list(range(len(df_marker_ranking)))

        return df_marker_ranking

    def get_marker_distribution_for_review(self, max_css: GUNC_RANKS, source: GuncRefDb, pct: float):

        # Obtain the markers that will actually be kept
        markers_at_pct_removed, markers_expected, domain_vote, markers_at_pct_removed_c = self.get_unq_markers_present_at_pct_removed(max_css, source, pct)
        markers_kept = frozenset(markers_at_pct_removed.keys())
        assert len(markers_kept) == len(markers_at_pct_removed)

        # We would like to find the last marker that was kept within that set
        contigs_markers_kept_is_on = {x['hit'].gene_id.rsplit('_')[0] for x in markers_at_pct_removed.values()}
        n_contigs_total = len(self.d_fna.keys())

        return len(contigs_markers_kept_is_on), n_contigs_total


    def split_fna_by_pct(self, max_css: GUNC_RANKS, source: GuncRefDb, pct: float):

        # Obtain the contig + marker ranking (and coordinates)
        df_marker_ranking = self.get_gunc_max_css_marker_ranking_and_coords(max_css, source)

        # Obtain the contig ranking (best to worst)
        contig_ranking = self.get_gunc_max_css_contig_ranking(max_css, source)
        contig_ranking = list(reversed([str(x) for x in contig_ranking['contig']]))

        # Load the fna file
        d_fna = self.d_fna

        total_nt = sum((len(x) for x in d_fna.values()))
        n_nt_to_keep = int(np.floor(total_nt * (pct / 100)))

        d_contig_to_coords_to_keep = dict()
        s_whole_contigs_to_keep = set()
        n_nt_kept = 0
        for contig in contig_ranking:
            cur_contig_len = len(d_fna[contig])

            # Stop processing if we've reached the desired length
            if n_nt_kept == n_nt_to_keep:
                break

            # We would take too many nucleotides, this contig needs to be split
            elif n_nt_kept + cur_contig_len > n_nt_to_keep:

                # Determine which part of this contig is the "better" half
                n_nt_to_keep_contig = n_nt_to_keep - n_nt_kept

                # Obtain the coordinates in each direction
                idx_fwd_from, idx_fwd_to = 0, n_nt_to_keep_contig
                idx_rev_from, idx_rev_to = cur_contig_len - n_nt_to_keep_contig, cur_contig_len

                # Subset the marker ranking dataframe to be only this contig
                df_subset = df_marker_ranking[df_marker_ranking['contig'] == contig]
                fwd_ranks = df_subset[
                    (df_subset['pos_from'] > idx_fwd_from) & (df_subset['pos_to'] <= idx_fwd_to)
                    ]['rank'].values
                rev_ranks = df_subset[
                    (df_subset['pos_from'] > idx_rev_from) & (df_subset['pos_to'] <= idx_rev_to)
                    ]['rank'].values

                # Take the direction that has the highest rank
                fwd_score = np.median(fwd_ranks) if len(fwd_ranks) > 0 else 0
                rev_score = np.median(rev_ranks) if len(rev_ranks) > 0 else 0

                # Calculate the index based on which direction is better
                if fwd_score >= rev_score:
                    d_contig_to_coords_to_keep[contig] = (
                        (0, n_nt_to_keep_contig), (n_nt_to_keep_contig, cur_contig_len)
                    )
                else:
                    d_contig_to_coords_to_keep[contig] = (
                        (cur_contig_len - n_nt_to_keep_contig, cur_contig_len),
                        (0, cur_contig_len - n_nt_to_keep_contig)
                    )

                # Stop processing as the remaining contigs will be discarded
                break

            # Otherwise, keep this contig and add more
            else:
                n_nt_kept += cur_contig_len
                s_whole_contigs_to_keep.add(contig)

        # Split the genome into fragment halves
        d_contig_to_seq_keep = dict()
        d_contig_to_seq_discard = dict()
        for contig in d_fna:
            if contig in d_contig_to_coords_to_keep:
                (idx_from_keep, idx_to_keep), (idx_from_disc, idx_to_disc) = d_contig_to_coords_to_keep[contig]
                d_contig_to_seq_keep[contig] = d_fna[contig][idx_from_keep:idx_to_keep]
                d_contig_to_seq_discard[contig] = d_fna[contig][idx_from_disc:idx_to_disc]
            elif contig in s_whole_contigs_to_keep:
                d_contig_to_seq_keep[contig] = d_fna[contig]
            else:
                d_contig_to_seq_discard[contig] = d_fna[contig]

        # Sanity check
        total_nt_kept = sum([len(x) for x in d_contig_to_seq_keep.values()])
        total_nt_disc = sum([len(x) for x in d_contig_to_seq_discard.values()])
        assert total_nt_kept + total_nt_disc == total_nt
        assert total_nt_kept == n_nt_to_keep

        return d_contig_to_seq_keep, d_contig_to_seq_discard

    def split_fna_by_pct_random(self, pct: float, seed: Optional[int] = None):

        # Obtain the random contig ranking (best to worst)
        contig_ranking = list(sorted(self.d_fna.keys()))
        random.Random(seed).shuffle(contig_ranking)

        # Load the fna file
        d_fna = self.d_fna

        total_nt = sum((len(x) for x in d_fna.values()))
        n_nt_to_keep = int(np.floor(total_nt * (pct / 100)))

        d_contig_to_coords_to_keep = dict()
        s_whole_contigs_to_keep = set()
        n_nt_kept = 0
        for contig in contig_ranking:
            cur_contig_len = len(d_fna[contig])

            # Stop processing if we've reached the desired length
            if n_nt_kept == n_nt_to_keep:
                break

            # We would take too many nucleotides, this contig needs to be split
            elif n_nt_kept + cur_contig_len > n_nt_to_keep:

                # Determine which part of this contig is the "better" half
                n_nt_to_keep_contig = n_nt_to_keep - n_nt_kept

                # Randomly take one half of this contig
                if random.Random(seed).random() < 0.5:
                    d_contig_to_coords_to_keep[contig] = (
                        (0, n_nt_to_keep_contig), (n_nt_to_keep_contig, cur_contig_len)
                    )
                else:
                    d_contig_to_coords_to_keep[contig] = (
                        (cur_contig_len - n_nt_to_keep_contig, cur_contig_len),
                        (0, cur_contig_len - n_nt_to_keep_contig)
                    )

                # Stop processing as the remaining contigs will be discarded
                break

            # Otherwise, keep this contig and add more
            else:
                n_nt_kept += cur_contig_len
                s_whole_contigs_to_keep.add(contig)

        # Split the genome into fragment halves
        d_contig_to_seq_keep = dict()
        d_contig_to_seq_discard = dict()
        for contig in d_fna:
            if contig in d_contig_to_coords_to_keep:
                (idx_from_keep, idx_to_keep), (idx_from_disc, idx_to_disc) = d_contig_to_coords_to_keep[contig]
                d_contig_to_seq_keep[contig] = d_fna[contig][idx_from_keep:idx_to_keep]
                d_contig_to_seq_discard[contig] = d_fna[contig][idx_from_disc:idx_to_disc]
            elif contig in s_whole_contigs_to_keep:
                d_contig_to_seq_keep[contig] = d_fna[contig]
            else:
                d_contig_to_seq_discard[contig] = d_fna[contig]

        # Sanity check
        total_nt_kept = sum([len(x) for x in d_contig_to_seq_keep.values()])
        total_nt_disc = sum([len(x) for x in d_contig_to_seq_discard.values()])
        assert total_nt_kept + total_nt_disc == total_nt
        assert total_nt_kept == n_nt_to_keep

        return d_contig_to_seq_keep, d_contig_to_seq_discard


def _test():
    GID = 'GCA_017985115.1'

    # df = pd.read_hdf('/srv/home/uqamussi/projects/gunc-chimeras/output/gunc/gunc_merged_max_css_level.h5')
    # row = df.loc[GID]

    g = Genome(GID)

    x = g.split_fna_by_pct(max_css='order', source=GuncRefDb.GTDB, pct=50)

    with open('/srv/home/uqamussi/tmp/fastani_keep.fna', 'w') as f:
        for contig, seq in x[0].items():
            f.write(f'>{contig}\n{seq}\n')
    with open('/srv/home/uqamussi/tmp/fastani_disc.fna' ,'w') as f:
        for contig, seq in x[1].items():
            f.write(f'>{contig}\n{seq}\n')


    return

    t = g.get_unq_markers_present_at_pct_removed(max_css='order', source=GuncRefDb.GTDB, pct=50)
    x, y = g.get_aligned_unq_domain_markers_from_markers(t, masked=True)

    loci = g.get_marker_loci()

    g.get_gunc_max_css_inferred_taxonomy('genus', source=GuncRefDb.GTDB)

    return

    g.get_marker_hits()
    loci = g.get_marker_loci()

    msa = g.get_domain_msa(omit_contigs=frozenset(t[0]))

    # y = g.get_random_contigs_where_removed_equals_x_pct_genome_removed(pct=50)
    #
    t_hits = g.get_marker_hits(omit_contigs=frozenset(t[0]))
    # y_hits = g.get_marker_hits(omit_contigs=frozenset(y[0]))
    return


if __name__ == '__main__':
    _test()
