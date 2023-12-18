import os
from typing import Set

import luigi
import pandas as pd
from chiton.fastani import fastani
from tqdm import tqdm
import multiprocessing as mp
import numpy as np
from workflow.config import DIR_OUT_FASTANI_CONTIG_SPLIT, DEBUG
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_sp_clusters import GtdbSpClustersR207
from workflow.fastani_contig_split.b_run_fastani import FastAniContigSplitRunFastAni
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.method.run_species_clustering_on_fastani_values import run_sp_clustering_on_fastani_values
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log
from workflow.util.paths import get_gid_root

import tempfile


class FastAniContigSplitZRunForSpecificCase(LuigiTask):
    target_pct = luigi.FloatParameter(default=50)

    query_gid = luigi.Parameter()

    @property
    def root_dir(self):
        return DIR_OUT_FASTANI_CONTIG_SPLIT

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    def output(self):
        return LocalTargetHdf5(
            os.path.join(
                self.root_dir, f'z_specific_case_{self.query_gid}_v6.h5'
            ))

    def run(self):
        log(f'FastAniContigSplitZRunForSpecificCase(p={self.target_pct}, gid={self.query_gid})', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()

        log('Loading MaxCSS scores')
        df_css = self.input()['max_css'].maybe_read_cached()

        # Set the genomes to run FastANI against
        if self.query_gid == 'GCA_900759445.1':
            ref_gids = {
                'g__Klebsiella': set(df_meta[(df_meta['genus'] == 'g__Klebsiella') & (df_meta['gtdb_representative'] == 't')].index),
                'g__Citrobacter': set(df_meta[(df_meta['genus'] == 'g__Citrobacter') & (df_meta['gtdb_representative'] == 't')].index)
            }
        else:
            raise Exception('TODO')

        css_row = df_css.loc[self.query_gid]

        # Write out the keep fna path
        with tempfile.TemporaryDirectory() as tmp_dir:
            path_keep_fna = os.path.join(tmp_dir, f'TEST_{self.query_gid}.fna')
            ref_gids = set.union(*ref_gids.values())

            # Write the keep fna
            write_genome_fna_worker(self.query_gid, self.target_pct, css_row['taxonomic_level'], css_row['source'], path_keep_fna)

            df = run_fastani_on_jobs(ref_gids, path_keep_fna, self.query_gid)

        df['query'] = self.query_gid
        # Save the result
        if not DEBUG:
            self.save_hdf(df)

        return


def write_genome_fna_worker(gid, pct_to_remove, gunc_tax_level, gunc_source, path_fna_keep):

    if gunc_source == 'gtdb':
        gunc_source = GuncRefDb.GTDB
    elif gunc_source == 'progenomes':
        gunc_source = GuncRefDb.PRO
    else:
        raise Exception(f'Unknown source: {gunc_source}')

    genome = Genome(gid)
    d_contig_to_seq_keep, d_contig_to_seq_disc = genome.split_fna_by_pct(
        max_css=gunc_tax_level,
        source=gunc_source,
        pct=pct_to_remove
    )

    # path_fna_disc = os.path.join(get_gid_root(gid, root=tmp_dir), f'{gid}_C.fna')
    # os.makedirs(os.path.dirname(path_fna_disc), exist_ok=True)

    # Write the modified FNA files to disk
    with open(path_fna_keep, 'w') as f:
        for contig, seq in d_contig_to_seq_keep.items():
            f.write(f'>{contig}\n{seq}\n')

    # with open(path_fna_disc, 'w') as f:
    #     for contig, seq in d_contig_to_seq_disc.items():
    #         f.write(f'>{contig}\n{seq}\n')

    return gid, path_fna_keep



def run_fastani_on_jobs(ref_gids: Set[str], path_fna_tmp, query_gid):
    ref_gid_paths = {x: os.path.join(get_gid_root(x), f'{x}.fna') for x in ref_gids}
    ref_gid_paths[f'TEST_{query_gid}'] = path_fna_tmp

    ani = fastani(query=list(ref_gid_paths.values()),
                  reference=list(ref_gid_paths.values()),
                  cpus=mp.cpu_count(),
                  single_execution=False,
                  bidirectional=True,
                  show_progress=True)

    # Prepare the ANI results for output
    keys = sorted(list(ref_gid_paths.keys()))
    arr_ani = np.zeros((len(keys), len(keys)), dtype=np.float64)
    arr_af = np.zeros((len(keys), len(keys)), dtype=np.float64)


    d_ani = ani.as_dict()

    for i, ref_gid in enumerate(keys):
        for j, qry_gid in enumerate(keys):

            q_key = ref_gid_paths[qry_gid]
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

            if q_key == r_key:
                ani = 100.0
                af = 1.0

            arr_ani[i, j] = ani
            arr_af[i, j] = af

    df = pd.DataFrame(arr_ani, index=keys, columns=keys)
    return df
