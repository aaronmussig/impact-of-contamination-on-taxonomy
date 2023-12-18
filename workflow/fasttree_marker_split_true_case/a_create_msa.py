import multiprocessing as mp
import os
import random
import tempfile
from collections import defaultdict, Counter

import luigi
import pandas as pd
from magna.util.disk import copy_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTTREE_MARKER_SPLIT_TRUE_CASE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_msa import GtdbMsaBacNonRepsR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetFasta, LocalTargetTsv
from workflow.model.taxonomy import TaxDomain
from workflow.util.log import log
from workflow.util.taxonomy import calculate_taxonomic_novelty

D_RANK = {'d': 'domain',
          'p': 'phylum',
          'c': 'class',
          'o': 'order',
          'f': 'family',
          'g': 'genus',
          's': 'species',
          'st': 'strain'}

RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
NOVELTY_TO_IDX = {'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6, 'strain': 7}


class FastTreeMarkerSplitTrueCaseCreateMsa(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)
    replicate_id = luigi.IntParameter()

    @property
    def root_dir(self):
        return os.path.join(
            DIR_OUT_FASTTREE_MARKER_SPLIT_TRUE_CASE,
            f'pct_{self.target_pct}',
            f'replicate_{self.replicate_id}'
        )

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
            'msa': GtdbMsaBacNonRepsR207(),
        }

    def output(self):
        return {
            'msa': LocalTargetFasta(os.path.join(self.root_dir, 'msa.faa')),
            'results': LocalTargetTsv(os.path.join(self.root_dir, 'results.tsv')),
            'replacements': LocalTargetTsv(os.path.join(self.root_dir, 'replacement_log.tsv'))
        }

    def load_true_msa(self, keep_gids):
        d_gid_to_msa = self.input()['msa'].maybe_read_cached()

        out = dict()
        for gid, seq in d_gid_to_msa.items():
            if gid[3:] in keep_gids:
                out[gid[3:]] = seq
        if len(out) != len(keep_gids):
            raise Exception(f'Expected {len(keep_gids)} gids, but only found {len(out)}')
        return out

    def run(self):
        log(f'FastTreeMarkerSplitTrueCaseCreateMsa(target_pct={self.target_pct}, r={self.replicate_id})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        bac_rep_gids = set(df_meta[(df_meta['gtdb_representative'] == 't') & (df_meta['domain'] == 'd__Bacteria')].index)
        bac_gids = set(df_meta[df_meta['domain'] == 'd__Bacteria'].index)
        d_gid_to_tax = df_meta['gtdb_taxonomy'].to_dict()

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        bac_gids_fail = bac_gids.intersection(set(df_max_css.index))

        # These are the genome ids that were used in the split tree (bac)
        original_split_gids = bac_gids_fail.union(bac_rep_gids)

        log(f'The original MSA would have contained: {len(original_split_gids):,}')
        log(f'Comprising of {len(bac_gids_fail):,} failed bacterial genomes')
        # bac_rep_gids_pass = rep_gids - bac_gids_fail

        log('Finding taxonomic novelty of failed genomes')
        d_og_gid_to_novelty, d_og_novelty_to_gid = calculate_taxonomic_novelty({x: d_gid_to_tax[x] for x in original_split_gids})

        # original_tax_novelty_counts = Counter([d_gid_to_novelty[x] for x in bac_gids_fail])
        # bac_gids_novelty_props = {k: 100 * v / n_gids_in_tree for k, v in bac_gids_novelty_counts.items()}

        # Populate this set with genomes to run on
        gids_to_run_on = bac_rep_gids - bac_gids_fail
        new_gids_added = set()
        gids_to_run_on_log = list()

        log('Creating a dictionary of taxonomy prefix to genomes')
        d_tax_prefix_to_gids = defaultdict(set)
        d_gids_to_run_on_split_tax = defaultdict(set)
        for gid, tax in tqdm(d_gid_to_tax.items()):
            tax_split = tax.split(';')
            for i in range(len(tax_split)):
                d_tax_prefix_to_gids[';'.join(tax_split[:i + 1])].add(gid)
                if gid in gids_to_run_on:
                    d_gids_to_run_on_split_tax[';'.join(tax_split[:i + 1])].add(gid)

        remaining_candidate_gids = bac_gids - bac_gids_fail - gids_to_run_on

        log('Finding replacement genomes for those which failed')
        for gid in tqdm(sorted(bac_gids_fail), smoothing=0.01):
            tax = d_gid_to_tax[gid]
            tax_split = tax.split(';')
            novelty = d_og_gid_to_novelty[gid]

            # Attempt to find a genome that satisfies the novelty criteria, that is closely related to the existing taxonomy
            start_idx = NOVELTY_TO_IDX[novelty]

            for i in reversed(range(len(tax_split))):
                tax_prefix = ';'.join(tax_split[:i + 1])

                # Are there any genomes from this tax string prefix that we can include?
                df_subset_gids = d_tax_prefix_to_gids[tax_prefix]
                gids_to_choose_from = list(remaining_candidate_gids.intersection(df_subset_gids))

                # We need to check what degree of taxonomic novelty each of the genomes we can choose from are
                # notably, you will have to check if adding it will alter the novelty of an existing genome

                # Will need to do a thorough check
                if novelty != 'strain':
                    new_gids_to_choose_from = set()
                    for candidate_gid in gids_to_choose_from:
                        candidate_tax = ';'.join(d_gid_to_tax[candidate_gid].split(';')[:start_idx+1])

                        n_existing_gids = len(d_gids_to_run_on_split_tax[candidate_tax])
                        if n_existing_gids == 0:
                            new_gids_to_choose_from.add(candidate_gid)
                    gids_to_choose_from = list(new_gids_to_choose_from)

                # Try and find a case where we can satisfy the novelty criterion
                if len(gids_to_choose_from) > 0:

                    replacement_gid = random.choice(gids_to_choose_from)
                    gids_to_run_on.add(replacement_gid)
                    new_gids_added.add(replacement_gid)
                    remaining_candidate_gids.remove(replacement_gid)

                    # Update the dictionary of taxonomic prefixes to genomes
                    tax_split_zz = d_gid_to_tax[replacement_gid].split(';')
                    for i_zz in range(len(tax_split_zz)):
                        d_gids_to_run_on_split_tax[';'.join(tax_split_zz[:i_zz + 1])].add(replacement_gid)

                    gids_to_run_on_log.append({
                        'gid': gid,
                        'gid_taxonomy': tax,
                        'replacement_gid': replacement_gid,
                        'replacement_gid_taxonomy': d_gid_to_tax[replacement_gid],
                        'tax_prefix': tax_split[i][0],
                        'novelty': novelty,
                    })
                    break
            else:
                # It's likely that there are no genomes that match the degree of taxonomic novelty
                # in this case, we would still like to make sure that we're including at least one genome
                # that is close to the taxonomy, run the same algorithm but don't limit on novelty
                for i in reversed(range(len(tax_split))):
                    tax_prefix = ';'.join(tax_split[:i + 1])

                    # Are there any genomes from this tax string prefix that we can include?
                    df_subset_gids = d_tax_prefix_to_gids[tax_prefix]
                    gids_to_choose_from = list(remaining_candidate_gids.intersection(df_subset_gids))

                    # Try and find a case where we can satisfy the novelty criterion
                    if len(gids_to_choose_from) > 0:
                        # df_subset = df_meta[df_meta.index.isin(gids_to_choose_from)]
                        # df_subset = df_subset.sort_values(by=['contig_count', 'checkm_completeness', 'checkm_contamination'], ascending=[True, False, True])
                        #
                        # replacement_gid = df_subset.index[0]

                        replacement_gid = random.choice(gids_to_choose_from)
                        gids_to_run_on.add(replacement_gid)
                        new_gids_added.add(replacement_gid)
                        remaining_candidate_gids.remove(replacement_gid)

                        # Update the dictionary of taxonomic prefixes to genomes
                        tax_split_zz = d_gid_to_tax[replacement_gid].split(';')
                        for i_zz in range(len(tax_split_zz)):
                            d_gids_to_run_on_split_tax[';'.join(tax_split_zz[:i_zz + 1])].add(replacement_gid)

                        gids_to_run_on_log.append({
                            'gid': gid,
                            'gid_taxonomy': tax,
                            'replacement_gid': replacement_gid,
                            'replacement_gid_taxonomy': d_gid_to_tax[replacement_gid],
                            'tax_prefix': tax_split[i][0],
                            'novelty': 'n/a',
                        })
                        break
                else:
                    raise Exception('Should not get here')

        log('Re-finding taxonomic novelty of failed genomes')
        d_new_gid_to_novelty, d_new_novelty_to_gid = calculate_taxonomic_novelty({x: d_gid_to_tax[x] for x in gids_to_run_on})

        for level in NOVELTY_TO_IDX:
            log(f'Number of genomes at {level} (old/new): {len(d_new_novelty_to_gid[level]):,} / {len(d_og_novelty_to_gid[level]):,}')

        df_replacement_log = pd.DataFrame(gids_to_run_on_log)
        log(f'Found {len(new_gids_added):,} new genomes to add, total={len(gids_to_run_on):,}')

        log('Creating queue')
        queue = [(x, self.target_pct) for x in sorted(new_gids_added)]
        if DEBUG:
            queue = queue[0:5]

        log('Processing queue')
        with mp.Pool(mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log('Collecting results')
        df_rows = list()
        d_qry_gid_to_seq = dict()
        for result in results:
            d_qry_gid_to_seq[result['gid']] = result['msa']
            d_qry_gid_to_seq[f"{result['gid']}_C"] = result['msa_c']
            cur_row = result.copy()
            del cur_row['msa']
            del cur_row['msa_c']
            df_rows.append(cur_row)
        df_results = pd.DataFrame(df_rows)

        log('Loading true MSA')
        d_true_msa = self.load_true_msa(gids_to_run_on)

        log(f'Processing the genomes that need to be run: {len(results):,}')
        lengths = set()
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp = os.path.join(tmp_dir, 'msa.faa')
            with open(path_tmp, 'w') as f:
                for gid, seq in sorted(d_true_msa.items()):
                    if gid in new_gids_added:
                        seq = d_qry_gid_to_seq[gid]
                        lengths.add(len(seq))
                        f.write(f'>TEST_{gid}\n{seq}\n')

                        seq = d_qry_gid_to_seq[f'{gid}_C']
                        lengths.add(len(seq))
                        f.write(f'>TEST_{gid}_C\n{seq}\n')

                    else:
                        lengths.add(len(seq))
                        f.write(f'>{gid}\n{seq}\n')

            if len(lengths) != 1:
                raise Exception(f'Expected all sequences to be the same length, but found {lengths}')

            if not DEBUG:
                log(f'Copying MSA ({get_file_size_fmt(path_tmp)})')
                copy_file(path_tmp, self.output()['msa'].path, checksum=True)

                path_csv_tmp = os.path.join(tmp_dir, 'results.tsv')
                df_results.to_csv(path_csv_tmp, sep='\t', index=False)
                log(f'Writing results ({get_file_size_fmt(path_csv_tmp)})')
                copy_file(path_csv_tmp, self.output()['results'].path, checksum=True)

                path_log_tmp = os.path.join(tmp_dir, 'results.tsv')
                df_replacement_log.to_csv(path_log_tmp, sep='\t', index=False)
                log(f'Writing replacement log ({get_file_size_fmt(path_log_tmp)}')
                copy_file(path_log_tmp, self.output()['replacements'].path, checksum=True)
        return


def worker(job):
    gid, target_pct = job
    results = {
        'gid': gid,
    }

    # Get the markers for this genome
    genome = Genome(gid)
    try:
        markers_at_pct_removed, markers_expected, domain, markers_at_pct_removed_c = genome.get_unq_random_markers_present_at_pct_removed(
            pct=target_pct)
    except KeyError:
        print(gid)
        raise

    # Sanity check
    if domain is not TaxDomain.BACTERIA:
        raise Exception(f'{gid} expected domain to be bacteria, but found {domain}')

    # Get the sequences of the markers at pct removed
    results['orig_markers'] = '|'.join(sorted(markers_expected))
    results['msa_markers'] = '|'.join(sorted(markers_at_pct_removed.keys()))
    results['msa_markers_c'] = '|'.join(sorted(markers_at_pct_removed_c.keys()))
    results['actual_pct'] = round(len(markers_at_pct_removed) / len(markers_expected) * 100, 4)

    # Create the MSA
    _, msa = genome.get_aligned_unq_domain_markers_from_markers(markers_at_pct_removed, masked=True)
    results['msa'] = msa

    # Create the MSA for the chimeric half
    _, msa_c = genome.get_aligned_unq_domain_markers_from_markers(markers_at_pct_removed_c, masked=True)
    results['msa_c'] = msa_c

    if len(msa) != len(msa_c) or len(msa) != 5036 or len(msa_c) != 5036:
        raise Exception('Expected all sequences to be 5036 bp')

    # Return the data
    return results


def get_taxonomic_novelty_v2(df_meta):
    # Get the number of taxa contained within each rank
    d_tax_under = defaultdict(lambda: 0)
    for taxonomy in df_meta['gtdb_taxonomy'].values:
        for rank in taxonomy.split(';'):
            d_tax_under[rank] += 1

    # Find the novelty
    n_total = len(df_meta)

    # Process the rows
    d_gid_to_tax_novelty = dict()
    for row in tqdm(df_meta.itertuples(), total=n_total):
        tax_str = row.gtdb_taxonomy
        taxonomy = tax_str.split(';')
        last_rank = 'st'
        for taxon in reversed(taxonomy):
            taxon_count = d_tax_under[taxon]
            if taxon_count > 1:
                break
            last_rank = taxon[0]
        d_gid_to_tax_novelty[row.Index] = D_RANK[last_rank]

    return d_gid_to_tax_novelty



