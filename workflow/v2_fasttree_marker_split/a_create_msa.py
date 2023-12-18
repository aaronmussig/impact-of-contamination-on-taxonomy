import multiprocessing as mp
import os
import tempfile

import luigi
import pandas as pd
from magna.util.disk import copy_file, get_file_size_fmt
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_V2_FASTTREE_MARKER_SPLIT
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_msa import GtdbMsaBacNonRepsR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetFasta, LocalTargetTsv
from workflow.util.log import log


class V2FastTreeMarkerSplitCreateMsa(LuigiTask):
    remove_pct = luigi.FloatParameter(default=50)

    @property
    def root_dir(self):
        return os.path.join(DIR_OUT_V2_FASTTREE_MARKER_SPLIT, f'remove_pct_{self.remove_pct}')

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
        log(f'V2FastTreeMarkerSplitCreateMsa(remove_pct={self.remove_pct})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        rep_gids = set(df_meta[(df_meta['gtdb_representative'] == 't') & (df_meta['domain'] == 'd__Bacteria')].index)

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        bac_rep_gids_fail = rep_gids.intersection(set(df_max_css.index))

        log('Loading true MSA')
        d_true_msa = self.load_true_msa(rep_gids)

        log('Creating queue')
        queue = list()
        queued_gids = set()
        for gid, row in df_max_css[df_max_css.index.isin(bac_rep_gids_fail)].iterrows():
            queued_gids.add(gid)
            queue.append((gid, row['source'], row['taxonomic_level'], 100 - self.remove_pct))
            if DEBUG and len(queue) > 5:
                break

        log('Processing queue')
        with mp.Pool(int(mp.cpu_count() * 0.8)) as pool:
            results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log('Collecting results')
        df_rows = list()
        d_qry_gid_to_seq = dict()
        for result in results:
            d_qry_gid_to_seq[result['gid']] = result['msa']
            cur_row = result.copy()
            del cur_row['msa']
            df_rows.append(cur_row)
        df_results = pd.DataFrame(df_rows)

        log(f'Processing the genomes that need to be run: {len(results):,}')
        lengths = set()
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_tmp = os.path.join(tmp_dir, 'msa.faa')
            with open(path_tmp, 'w') as f:
                for gid, seq in sorted(d_true_msa.items()):
                    if gid in queued_gids:
                        seq = d_qry_gid_to_seq[gid]
                        lengths.add(len(seq))
                        f.write(f'>TEST_{gid}\n{seq}\n')
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
        return


def worker(job):
    gid, source, max_css, target_pct = job
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

    # Get the sequences of the markers at pct removed
    markers_at_pct_removed, markers_expected, domain_vote, markers_at_pct_removed_c = genome.get_unq_markers_present_at_pct_removed(
        max_css=max_css, source=source, pct=target_pct
    )
    results['orig_markers'] = '|'.join(sorted(markers_expected))
    results['msa_markers'] = '|'.join(sorted(markers_at_pct_removed.keys()))
    results['actual_pct_removed'] = round(100 - (len(markers_at_pct_removed) / len(markers_expected) * 100), 4)

    # Create the MSA
    _, msa = genome.get_aligned_unq_domain_markers_from_markers(markers_at_pct_removed, masked=True)
    results['msa'] = msa

    # Return the data
    return results
