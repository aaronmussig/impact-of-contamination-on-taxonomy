import multiprocessing as mp
import os

import luigi
import pandas as pd
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_REVIEW
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class ReviewMarkerDistribution(LuigiTask):
    target_pct = luigi.FloatParameter()

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTarget(
            os.path.join(
                DIR_OUT_REVIEW,
                f'marker_distribution__p{self.target_pct}.tsv'
            ))

    def run(self):
        log(f'ReviewMarkerDistribution(p={self.target_pct})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        rep_gids = set(df_meta[(df_meta['gtdb_representative'] == 't') & (df_meta['domain'] == 'd__Bacteria')].index)

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        bac_rep_gids_fail = rep_gids.intersection(set(df_max_css.index))

        log('Creating queue')
        queue = list()
        queued_gids = set()
        for gid, row in df_max_css[df_max_css.index.isin(bac_rep_gids_fail)].iterrows():
            queued_gids.add(gid)
            queue.append((gid, row['source'], row['taxonomic_level'], 100 - self.target_pct))
            if DEBUG and len(queue) > 5:
                break

        log('Processing queue')
        if DEBUG:
            results = list(worker(x) for x in queue)
        else:
            with mp.Pool(int(mp.cpu_count() * 0.8)) as pool:
                results = list(tqdm(pool.imap_unordered(worker, queue), total=len(queue)))

        log(f'Writing to: {self.output().path}')
        df = pd.DataFrame(results)
        df.to_csv(self.output().path, index=False, sep='\t')

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

    # Get the marker ranking
    n_contigs_used, n_contigs_total = genome.get_marker_distribution_for_review(max_css, source, pct=50)
    results['n_contigs_used'] = n_contigs_used
    results['n_contigs_total'] = n_contigs_total

    # Return the data
    return results
