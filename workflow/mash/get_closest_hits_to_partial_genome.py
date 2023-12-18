import multiprocessing as mp
import os
import tempfile
from collections import defaultdict
from typing import Set

import luigi
import pandas as pd
from chiton.fastani import fastani
from magna.util.disk import copy_file, move_file
from tqdm import tqdm

from workflow.config import DIR_OUT_FASTANI_CONGRUENCE, DEBUG, DIR_OUT_BATCH, DIR_CACHE
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fastani_congruence.a_run_mash_on_failed import FastAniCongruenceRunMash
from workflow.genome.closest_100_genomes_to_representative import Closest100GenomesToRepresentative
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.gunc_helper.genome_pct_congruence_contigs_removed import GenomePctCongruenceContigsRemoved
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.collection import iter_batches
from workflow.util.log import log
from workflow.util.paths import get_gid_root
from workflow.util.rq import rq_wait_for_queue_empty, rq_and_wait


class MashGetClosestHitsToPartialGenome(LuigiTask):
    target_pct = luigi.FloatParameter()
    congruence = luigi.FloatParameter()

    mash_k = luigi.IntParameter(default=16)
    mash_s = luigi.IntParameter(default=5000)

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'max_css': AggregateMaxCssLevelMerged(),
        }

    # def output(self):
    #     return LocalTargetHdf5(
    #         os.path.join(
    #             DIR_OUT_FASTANI_CONGRUENCE,
    #             f'ani_k{self.mash_k}_s{self.mash_s}__c{self.congruence}_p{self.target_pct}.h5'
    #         ))

    def run(self):
        log(f'MashGetClosestHitsToPartialGenome(p={self.target_pct}, c={self.congruence})', title=True)
        # self.make_output_dirs()

        gid = 'GCF_013361095.1'


        ref_genome = Genome(gid)

        # Get the congruence value of the reference genome
        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        ref_css = df_max_css.loc[gid, 'taxonomic_level']
        ref_source = df_max_css.loc[gid, 'source']
        ref_source = GuncRefDb(ref_source)

        df_contig_ranking = ref_genome.get_gunc_max_css_contig_ranking(ref_css, ref_source)
        contig_order = list(df_contig_ranking['contig'].values)

        with open('/tmp/genome.fna', 'w') as f:
            for contig in contig_order:
                if contig == 'NZ_JAATWA010000006.1':
                    break
                seq = ref_genome.d_fna[contig]
                f.write(f'>{contig}\n{seq}\n')

        return


