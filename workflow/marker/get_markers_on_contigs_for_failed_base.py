import multiprocessing as mp
import os
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_MARKER, R207_MARKERS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid_include_mul
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.paths import get_gid_root


class GetMarkersOnContigsForFailedBase(LuigiTask):
    """Store the gid / contig / number of R207 markers on the contig"""

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'markers_on_contigs_for_failed_base.h5'))

    def run(self):
        log('Getting markers on contigs for GUNC failed genomes (base)', title=True)
        self.make_output_dirs()

        log('Loading Max CSS')
        df_css = self.input()['max_css'].read() if not DEBUG else self.input()['max_css'].read_cached()

        log('Creating queue')
        queue = sorted(set(df_css.index))

        log('Processing queue')
        if DEBUG:
            results = [run_on_gid(x) for x in queue[:5]]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = list(tqdm(pool.imap_unordered(run_on_gid, queue), total=len(queue)))

        log('Creating dataframe')
        rows = list()
        for result in results:
            rows.extend(result)
        df = pd.DataFrame(rows)

        log('Saving dataframe')
        if not DEBUG:
            self.save_hdf(df, index=['gid', 'contig'])


def run_on_gid(gid):
    # Read the FASTA file
    gid_root = get_gid_root(gid)

    # Load the top hit files
    path_faa = os.path.join(gid_root, 'prodigal', f'{gid}.faa')
    d_faa = dict(read_fasta(path_faa))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # Iterate over each percentage to determine which makers are kept
    _, marker_list = get_marker_hits_for_gid_include_mul(d_faa, pfam_th, tigr_th)

    d_contig_to_marker = defaultdict(lambda: defaultdict(lambda: 0))
    for marker in marker_list:
        cur_hit = marker['hit']
        contig_id = cur_hit.gene_id[0:cur_hit.gene_id.rindex('_')]
        d_contig_to_marker[contig_id][cur_hit.hmm_id] += 1

    out = list()
    for contig_id, d_marker_to_cnt in d_contig_to_marker.items():
        cur_row = {
            'gid': gid,
            'contig': contig_id
        }
        for marker in R207_MARKERS:
            cur_row[marker] = d_marker_to_cnt[marker]
        out.append(cur_row)
    return out
