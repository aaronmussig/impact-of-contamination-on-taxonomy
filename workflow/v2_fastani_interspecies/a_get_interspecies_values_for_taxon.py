import os
import tempfile
from typing import Set

import luigi

from workflow.config import DIR_OUT_V2_FASTANI_INTERSP, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
import multiprocessing as mp

import luigi
import pandas as pd
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.util.paths import get_gid_root


class V2FastAniInterSpeciesGetValuesForTaxon(LuigiTask):
    taxon = luigi.Parameter()

    @property
    def root_dir(self):
        return DIR_OUT_V2_FASTANI_INTERSP

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(self.root_dir, f'{self.taxon}.h5'))

    def run(self):
        log(f'V2FastAniInterSpeciesGetValuesForTaxon(taxon={self.taxon})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        gids_in_taxon = set()
        for rank in ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'):
            gids_in_taxon.update(set(df_meta[df_meta[rank] == self.taxon].index))
        rep_gids = set(df_meta[df_meta['gtdb_representative'] == 't'].index)
        rep_gids_in_taxon = rep_gids.intersection(gids_in_taxon)
        non_rep_gids_in_taxon = gids_in_taxon.difference(rep_gids_in_taxon)

        log(f'Found {len(rep_gids_in_taxon):,} representative genomes in taxon')
        log(f'Found {len(non_rep_gids_in_taxon):,} non-representative genomes in taxon')

        df = run_fastani_on_jobs(rep_gids_in_taxon, non_rep_gids_in_taxon)

        if not DEBUG:
            self.save_hdf(df)

        return


def run_fastani_on_jobs(ref_gids: Set[str], qry_gids: Set[str]):
    ref_gid_paths = {x: os.path.join(get_gid_root(x), f'{x}.fna') for x in ref_gids}
    qry_gid_paths = {x: os.path.join(get_gid_root(x), f'{x}.fna') for x in qry_gids}

    out = list()

    ani = fastani(query=list(qry_gid_paths.values()),
                  reference=list(ref_gid_paths.values()),
                  cpus=mp.cpu_count(),
                  single_execution=False,
                  bidirectional=True,
                  show_progress=True)

    # Prepare the ANI results for output
    d_ani = ani.as_dict()
    for qry_gid in qry_gid_paths.keys():
        q_key = qry_gid_paths[qry_gid]
        for ref_gid in ref_gid_paths.keys():
            r_key = ref_gid_paths[ref_gid]
            qvr = d_ani[q_key][r_key]
            rvq = d_ani[r_key][q_key]

            if qvr is not None and rvq is not None:
                ani = max(qvr.ani, rvq.ani)
                af = max(qvr.align_frac, rvq.align_frac)
            elif qvr is not None and rvq is None:
                ani = qvr.ani
                af = qvr.align_frac
            elif qvr is None and rvq is not None:
                ani = rvq.ani
                af = rvq.align_frac
            else:
                ani = 0
                af = 0

            out.append({
                'query': qry_gid,
                'ref': ref_gid,
                'ani': ani,
                'af': af,
            })

    df = pd.DataFrame(out)
    return df
