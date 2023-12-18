import multiprocessing as mp
import os
import random
import tempfile

import luigi
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM
from workflow.exception import GenomeNoMarkersFound
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_msa import GtdbMsaBacR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.luigi import LuigiTask, LocalTargetFasta
from workflow.model.taxonomy import TaxDomain
from workflow.util.log import log


class FastTreeMarkerAnalysisBatchesRandom(LuigiTask):
    target_pct = luigi.FloatParameter()
    batch_id = luigi.IntParameter()

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
            'msa': GtdbMsaBacR207(),
        }

    def output(self):
        return LocalTargetFasta(
            os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, f'{self.batch_id}_{self.target_pct}', 'msa.faa'))

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
        log(f'Creating batches for a random selection of reps (batch={self.batch_id} @ {self.target_pct}%)', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        rep_gids = set(df_meta[(df_meta['gtdb_representative'] == 't') & (df_meta['domain'] == 'd__Bacteria')].index)

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        fail_gids = set(df_max_css.index)
        n_failed_reps = len(fail_gids & rep_gids)

        log('Loading true MSA')
        d_true_msa = self.load_true_msa(rep_gids)

        log('Determining random gids to use')
        gids_to_use = set(random.sample(rep_gids, n_failed_reps))
        assert len(gids_to_use) == n_failed_reps

        log('Creating queue')
        queue = list()
        for gid in sorted(gids_to_use):
            queue.append((gid, self.target_pct))

            if DEBUG and len(queue) > 5:
                break

        log(f'Processing {len(queue):,} gids')
        with mp.Pool(processes=mp.cpu_count() if not DEBUG else 1) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log('Processing results to check if any results do not need to be run')
        gids_changed_domain = set()
        gids_same_msa = set()
        gids_to_run = dict()
        gids_no_markers = set()
        for result in results:
            gid = result['gid']
            expected_domain = df_meta.loc[gid]['domain']

            if result['domain'] is None:
                gids_no_markers.add(gid)
                log(f'No markers for gid: {gid}')
                continue

            # Case 1: Genome changed domains (due to loss of marker genes, this will be wrong!)
            if result['domain'] != 'd__Bacteria' or expected_domain != result['domain']:
                gids_changed_domain.add(gid)
                log(f'Expected domain {expected_domain} != {result["domain"]} for {gid}')
                continue

            # Case 2: MSA has not changed (this is correct)
            if result['msa'] == d_true_msa[gid]:
                gids_same_msa.add(gid)
                log(f'MSA for {gid} has not changed')
                continue

            # Case 3: No markers were found (will be incorrect)
            if set(result['msa']) == {'-'}:
                raise Exception(f'No markers found for {gid}')

            # Otherwise, this needs to be run
            gids_to_run[gid] = result

        log(f'Processing the genomes that need to be run: {len(gids_to_run):,}')
        lengths = set()
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp = os.path.join(tmp_dir, 'msa.faa')
            with open(path_tmp, 'w') as f:
                for gid, seq in sorted(d_true_msa.items()):
                    if gid in gids_to_run:
                        seq = gids_to_run[gid]['msa']
                        gid = f'TEST_{gid}'
                    lengths.add(len(seq))
                    f.write(f'>{gid}\n{seq}\n')

            if len(lengths) != 1:
                raise Exception(f'Expected all sequences to be the same length, but found {lengths}')

            root_dir = os.path.join(DIR_OUT_FASTTREE_MARKER_ANALYSIS_RANDOM, f'{self.batch_id}_{self.target_pct}')
            with open(os.path.join(root_dir, 'gids_changed_domain.txt'), 'w') as f:
                f.write('\n'.join(sorted(gids_changed_domain)))
                f.write('\n')

            with open(os.path.join(root_dir, 'gids_same_msa.txt'), 'w') as f:
                f.write('\n'.join(sorted(gids_same_msa)))
                f.write('\n')

            with open(os.path.join(root_dir, 'gids_no_markers.txt'), 'w') as f:
                f.write('\n'.join(sorted(gids_no_markers)))
                f.write('\n')

            with open(os.path.join(root_dir, 'gids_to_run_on.txt'), 'w') as f:
                f.write('\n'.join(sorted(gids_to_run.keys())))
                f.write('\n')

            with open(os.path.join(root_dir, 'gids_contigs_removed.txt'), 'w') as f:
                for result in sorted(results, key=lambda x: x['gid']):
                    contigs_removed_str = ';'.join(sorted(result['contigs']))
                    f.write(f'{result["gid"]}\t{result["pct_removed"]}\t{contigs_removed_str}\n')

            log('Copying file')
            if not DEBUG:
                copy_file(path_tmp, self.output().path, checksum=True)
        return


def worker(job):
    gid, target_pct = job
    results = {
        'gid': gid
    }

    # Randomly remove the contigs
    genome = Genome(gid)
    contigs, pct_removed = genome.get_random_contigs_where_removed_equals_x_pct_genome_removed(target_pct)
    contigs = frozenset(contigs)
    results['contigs'] = contigs
    results['pct_removed'] = pct_removed

    # Calculate the new domain of the genome using the markers that remain
    try:
        d_domain_results = genome.get_domain_from_markers(omit_contigs=contigs)
    except GenomeNoMarkersFound:
        results['domain'] = None
        return results
    domain = d_domain_results['domain']
    if domain is not TaxDomain.BACTERIA:
        results['domain'] = 'd__Archaea'
        return results
    results['domain'] = 'd__Bacteria'

    # Otherwise, get the MSA
    msa = genome.get_domain_msa(masked=True, omit_contigs=contigs)
    results['msa'] = msa
    return results
