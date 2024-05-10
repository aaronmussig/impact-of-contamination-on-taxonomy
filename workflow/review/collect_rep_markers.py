import multiprocessing as mp
import os

import luigi
from luigi import LocalTarget
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_REVIEW, \
    DIR_OUT_REVIEW_CHECKM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask
from workflow.util.log import log


class ReviewCollectRepMarkers(LuigiTask):
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
                f'marker_distribution__p{self.target_pct}.123'
            ))

    def run(self):
        log(f'ReviewMarkerDistribution(p={self.target_pct})', title=True)
        self.make_output_dirs()
        os.makedirs(DIR_OUT_REVIEW_CHECKM, exist_ok=True)

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

    out_path = os.path.join(DIR_OUT_REVIEW_CHECKM, f'{gid}.faa')
    with open(out_path, 'w') as f:
        for marker, d_hit in markers_at_pct_removed.items():
            cur_hit = d_hit['hit']
            cur_gene = cur_hit.gene_id
            cur_seq = d_hit['seq']
            f.write(f'>{cur_gene}_{marker}\n{cur_seq}\n')

    # checkm lineage_wf -t 90 --genes -x faa /srv/home/uqamussi/projects/gunc-chimeras/output/review/checkm_input /srv/home/uqamussi/projects/gunc-chimeras/output/review/checkm_output

    # Return the data
    return results
