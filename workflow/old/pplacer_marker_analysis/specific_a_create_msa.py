import multiprocessing as mp
import os
import tempfile
from typing import Dict, Set

import luigi
import pandas as pd
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_PPLACER_MARKER_ANALYSIS
from workflow.exception import GenomeNoMarkersFound
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_msa import GtdbMsaBacNonRepsR207, GtdbMsaArcNonRepsR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetFasta, LocalTargetHdf5
from workflow.model.taxonomy import TaxDomain
from workflow.util.log import log


class PplacerMarkerAnalysisSpecificCreateMsa(LuigiTask):
    pplacer_domain = luigi.ChoiceParameter(choices=['arc', 'bac'])
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()
    gids = luigi.ListParameter()

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
            'msa_bac': GtdbMsaBacNonRepsR207(),
            'msa_arc': GtdbMsaArcNonRepsR207(),
        }

    @property
    def root_dir(self) -> str:
        return os.path.join(
            DIR_OUT_PPLACER_MARKER_ANALYSIS,
            f'specific_{self.pplacer_domain}_p{self.target_pct}_c{self.congruence}__{"_".join(sorted(self.gids))}'
        )

    def output(self):
        return {
            'msa': LocalTargetFasta(os.path.join(self.root_dir, 'msa.fasta')),
            'log': LocalTargetHdf5(os.path.join(self.root_dir, 'log.h5')),
        }

    def load_true_msa(self, keep_gids: Set[str]) -> Dict[str, str]:
        d_gid_to_msa = {
            **self.input()['msa_bac'].maybe_read_cached(),
            **self.input()['msa_arc'].maybe_read_cached()
        }
        out = dict()
        for gid, seq in d_gid_to_msa.items():
            if gid[3:] in keep_gids:
                out[gid[3:]] = seq
        if len(out) != len(keep_gids):
            raise Exception(f'Expected {len(keep_gids)} gids, but only found {len(out)}')
        return out

    def run(self):
        log(f'PplacerMarkerAnalysisSpecificCreateMsa('
            f'congruence={self.congruence}, pct={self.target_pct}, '
            f'domain={self.pplacer_domain}, gids={self.gids})', title=True)
        self.make_output_dirs()

        # Convert the domain to d__Archaea or d__Bacteria
        domain_str = _domain_to_gtdb_str(self.pplacer_domain)

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        df_meta_subset = df_meta[df_meta.index.isin(self.gids)]
        assert (len(df_meta_subset) == len(self.gids))
        log(f'Matched all {len(df_meta_subset):,} input genomes to the metadata')

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()

        log('Subset dataframe to only those fail non rep gids')
        df_max_css_subset = df_max_css[df_max_css.index.isin(self.gids)]
        log(f'DF CSS Length: {len(df_max_css_subset):,}')

        log('Loading true MSA for gids to run on')
        d_true_msa = self.load_true_msa(set(self.gids))

        log('Creating queue')
        queue = create_queue_gunc(df_max_css_subset, self.congruence, self.target_pct)

        log('Processing queue')
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(work_queue_gunc, queue), total=len(queue), smoothing=0.01))

        log('Processing results to check if any results do not need to be run')
        rows = list()
        gids_to_run = dict()
        for result in results:
            gid = result['gid']
            expected_domain = domain_str
            cur_row = {
                'gid': gid,
                'domain': result['domain'],
                'result': 'N/A',
                'pct_removed': result['pct_removed'],
                'contigs': ';'.join(result['contigs']),
                'original_markers': ';'.join(result['original_markers']),
                'new_markers': ';'.join(result['new_markers']),
            }

            if result['domain'] is None:
                cur_row['result'] = 'no_markers'
                cur_row['domain'] = 'N/A'
                rows.append(cur_row)
                continue

            # Case 1: Genome changed domains (due to loss of marker genes, this will be wrong!)
            if expected_domain != result['domain']:
                cur_row['result'] = 'changed_domain'
                rows.append(cur_row)
                continue

            # Case 2: MSA has not changed (this is correct)
            if result['msa'] == d_true_msa[gid]:
                cur_row['result'] = 'same_msa'
                rows.append(cur_row)
                continue

            # Case 3: No markers were found (will be incorrect)
            if set(result['msa']) == {'-'}:
                raise Exception(f'No markers found for {gid}')

            # Otherwise, this needs to be run
            gids_to_run[gid] = result
            cur_row['result'] = 'run'
            rows.append(cur_row)

        # Write the log file
        assert (len(rows) == len(queue))
        df_log = pd.DataFrame(rows)
        if not DEBUG:
            self.save_hdf(df_log, index='gid', path=self.output()['log'].path)

        # Write the MSA
        log(f'Creating MSA for: {len(gids_to_run):,}')
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp_msa = os.path.join(tmp_dir, 'msa.faa')

            with open(path_tmp_msa, 'w') as f_msa:
                for gid, cur_result in sorted(gids_to_run.items()):
                    f_msa.write(f'>TEST_{gid}\n{cur_result["msa"]}\n')

            path_msa = self.output()['msa'].path

            log('Copying file')
            if not DEBUG:
                copy_file(path_tmp_msa, path_msa, checksum=True)
        return


def worker(job):
    gid, source, max_css, congruence, target_pct = job
    results = {
        'gid': gid,
        'source': source,
        'max_css': max_css
    }

    if source == 'gtdb':
        source = GuncRefDb.GTDB
    elif source == 'progenomes':
        source = GuncRefDb.PRO
    else:
        raise ValueError(f'Unknown source: {source}')

    genome = Genome(gid)
    contigs, pct_removed = genome.get_gunc_contigs_where_removed_equals_x_pct_genome_removed(
        pct=target_pct,
        max_congruence=congruence,
        max_css=max_css,
        source=source
    )
    results['pct_removed'] = pct_removed
    contigs = frozenset(contigs)

    # Calculate the new domain of the genome using the markers that remain
    try:
        d_domain_results = genome.get_domain_from_markers(omit_contigs=contigs)
    except GenomeNoMarkersFound:
        results['domain'] = None
        return results

    domain = d_domain_results['domain']
    if domain is TaxDomain.BACTERIA:
        results['domain'] = 'd__Bacteria'
    elif domain is TaxDomain.ARCHAEA:
        results['domain'] = 'd__Archaea'
    else:
        raise ValueError(f'Unknown domain: {domain}')

    # Otherwise, get the MSA
    msa = genome.get_domain_msa(masked=True, omit_contigs=contigs)
    results['msa'] = msa
    return results


def _domain_to_gtdb_str(domain):
    if domain == 'arc':
        return 'd__Archaea'
    elif domain == 'bac':
        return 'd__Bacteria'
    else:
        raise ValueError(f'Invalid domain: {domain}')


def create_queue_random(gids_to_run_on: Set[str], target_pct: float):
    queue = list()
    for gid in gids_to_run_on:
        queue.append((gid, target_pct))
        if DEBUG and len(queue) > 5:
            break
    return queue


def create_queue_gunc(df_max_css: pd.DataFrame, congruence: float, target_pct: float):
    queue = list()
    for gid, row in df_max_css.iterrows():
        queue.append((gid, row['source'], row['taxonomic_level'], congruence, target_pct))
    return queue


def work_queue_random(job):
    gid, target_pct = job
    results = {
        'gid': gid,
    }

    genome = Genome(gid)
    contigs, pct_removed = genome.get_random_contigs_where_removed_equals_x_pct_genome_removed(
        pct=target_pct,
    )
    results['pct_removed'] = pct_removed
    contigs = frozenset(contigs)
    results['contigs'] = contigs

    # Calculate the original markers
    genome_original_markers = genome.get_marker_hits()
    results['original_markers'] = frozenset({**genome_original_markers['unq'], **genome_original_markers['muq']}.keys())

    genome_new_markers = genome.get_marker_hits(omit_contigs=contigs)
    results['new_markers'] = frozenset({**genome_new_markers['unq'], **genome_new_markers['muq']}.keys())

    # Calculate the new domain of the genome using the markers that remain
    try:
        d_domain_results = genome.get_domain_from_markers(omit_contigs=contigs)
    except GenomeNoMarkersFound:
        results['domain'] = None
        return results

    domain = d_domain_results['domain']
    if domain is TaxDomain.BACTERIA:
        results['domain'] = 'd__Bacteria'
    elif domain is TaxDomain.ARCHAEA:
        results['domain'] = 'd__Archaea'
    else:
        raise ValueError(f'Unknown domain: {domain}')

    # Otherwise, get the MSA
    msa = genome.get_domain_msa(masked=True, omit_contigs=contigs)
    results['msa'] = msa
    return results


def work_queue_gunc(job):
    gid, source, max_css, congruence, target_pct = job
    results = {
        'gid': gid,
        'source': source,
        'max_css': max_css
    }

    if source == 'gtdb':
        source = GuncRefDb.GTDB
    elif source == 'progenomes':
        source = GuncRefDb.PRO
    else:
        raise ValueError(f'Unknown source: {source}')

    genome = Genome(gid)
    contigs, pct_removed = genome.get_gunc_contigs_where_removed_equals_x_pct_genome_removed(
        pct=target_pct,
        max_congruence=congruence,
        max_css=max_css,
        source=source
    )
    results['pct_removed'] = pct_removed
    contigs = frozenset(contigs)
    results['contigs'] = contigs

    # Calculate the original markers
    genome_original_markers = genome.get_marker_hits()
    results['original_markers'] = frozenset({**genome_original_markers['unq'], **genome_original_markers['muq']}.keys())

    genome_new_markers = genome.get_marker_hits(omit_contigs=contigs)
    results['new_markers'] = frozenset({**genome_new_markers['unq'], **genome_new_markers['muq']}.keys())

    # Calculate the new domain of the genome using the markers that remain
    try:
        d_domain_results = genome.get_domain_from_markers(omit_contigs=contigs)
    except GenomeNoMarkersFound:
        results['domain'] = None
        return results

    domain = d_domain_results['domain']
    if domain is TaxDomain.BACTERIA:
        results['domain'] = 'd__Bacteria'
    elif domain is TaxDomain.ARCHAEA:
        results['domain'] = 'd__Archaea'
    else:
        raise ValueError(f'Unknown domain: {domain}')

    # Otherwise, get the MSA
    msa = genome.get_domain_msa(masked=True, omit_contigs=contigs)
    results['msa'] = msa
    return results
