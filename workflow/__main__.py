import os
import shutil
import socket
import time

import luigi as lg
import typer
from redis import Redis
from rq import Queue
from rq.command import send_shutdown_command
from rq.worker import Worker

from workflow.blast.a_blast_genome_against_207_db import BlastGenomeAgainstR207Db
from workflow.blast.b_get_taxstring_against_hits import BlastGetTaxStringAgainstHits
from workflow.circos_plot.create_circos_plot import CreateCircosPlot
from workflow.config import REDIS_HOST, REDIS_PASS, DIR_CACHE
from workflow.external.gtdb_r207_blastdb import GtdbR207BlastDb
from workflow.fastani_contig_split.a_run_mash import FastAniContigSplitRunMash
from workflow.fastani_contig_split.b_run_fastani import FastAniContigSplitRunFastAni
from workflow.fastani_contig_split.c_analyse_fastani import FastAniContigSplitAnalyseFastAni
from workflow.fastani_contig_split.d_report_results import FastAniContigSplitReportResultsFastAni
from workflow.fastani_contig_split.d_report_results_relaxed import FastAniContigSplitReportResultsFastAniRelaxed
from workflow.fastani_contig_split.random_a_run_mash import FastAniContigSplitRandomRunMash
from workflow.fastani_contig_split.random_b_run_fastani import FastAniContigSplitRandomRunFastAni
from workflow.fastani_contig_split.random_c_analyse_fastani import FastAniContigSplitRandomAnalyseFastAni
from workflow.fastani_contig_split.random_d_report_results import FastAniContigSplitRandomReportResultsFastAni
from workflow.fastani_contig_split.z_run_for_specific_case import FastAniContigSplitZRunForSpecificCase
from workflow.fastani_random.e_analyse_results import FastAniRandomAnalyseResults
from workflow.fasttree_marker_split.a_create_msa import FastTreeMarkerSplitCreateMsa
from workflow.fasttree_marker_split.b_run_fasttree import FastTreeMarkerSplitRunFastTree
from workflow.fasttree_marker_split.c_root_fasttree import FastTreeMarkerSplitRootFastTree
from workflow.fasttree_marker_split.e_decorate_fasttree import FastTreeMarkerSplitDecorateFastTree
from workflow.fasttree_marker_split.f_a_compare_patristic_distance import FastTreeMarkerSplitComparePatristicDistance
from workflow.fasttree_marker_split.f_analyse_decorated import FastTreeMarkerSplitDecorateAnalyseFastTree
from workflow.fasttree_marker_split.f_b_calculate_num_jumps_between_halves import FastTreeMarkerSplitJumpsBetweenHalves
from workflow.fasttree_marker_split.random_a_create_msa import FastTreeMarkerSplitRandomCreateMsa
from workflow.fasttree_marker_split.random_d_decorate_fasttree import FastTreeMarkerSplitRandomDecorateFastTree
from workflow.fasttree_marker_split_true_case.a_create_msa import FastTreeMarkerSplitTrueCaseCreateMsa
from workflow.fasttree_marker_split_true_case.c_calculate_patristic_distance import \
    FastTreeMarkerSplitTrueCaseCalculatePatristicDistance
from workflow.fasttree_randomly_remove_markers.a_create_msa import FastTreeRandomlyRemoveMarkersCreateMsa
from workflow.fasttree_randomly_remove_markers.c_root_fasttree import FastTreeRandomlyRemoveMarkersRootFastTree
from workflow.fasttree_randomly_remove_markers.d_decorate_fasttree import FastTreeRandomlyRemoveMarkersDecorateFastTree
from workflow.fasttree_randomly_remove_markers.e_analyse_decorated import FastTreeRandomlyRemoveMarkersAnalyseDecorated
from workflow.final.generate_master_tsv import FinalGenerateMasterTsv
from workflow.final.v2_generate_master_tsv import V2FinalGenerateMasterTsv
from workflow.marker_trimming.gunc_f_analyse_decorated import MarkerTrimmingGuncDecorateAnalyseFastTree
from workflow.v2_fastani_interspecies.a_get_interspecies_values_for_taxon import V2FastAniInterSpeciesGetValuesForTaxon
from workflow.v2_fastani_interspecies.b_get_interspecies_for_failed_to_genus import \
    V2FastAniInterSpeciesGetValuesForFailedToGenus
from workflow.v2_fastani_rep_to_closest_rep.a_create_jobs import V2FastAniRepToClosestRepCreateJobs
from workflow.v2_fasttree_marker_split.a_create_msa import V2FastTreeMarkerSplitCreateMsa
from workflow.v2_fasttree_marker_split.c_root_fasttree import V2FastTreeMarkerSplitRootFastTree
from workflow.v2_fasttree_marker_split.d_decorate_fasttree import V2FastTreeMarkerSplitDecorateFastTree
from workflow.v2_fasttree_marker_split.e_analyse_decorated import V2FastTreeMarkerSplitAnalyseDecoratedFastTree
from workflow.v2_pplacer_marker_split.a_create_msa import V2PplacerMarkerSplitCreateMsa
from workflow.v2_pplacer_marker_split.b_run_pplacer import V2PplacerMarkerSplitRunPplacer
from workflow.v2_pplacer_marker_split.c_get_taxonomy_from_tree import V2PplacerMarkerSplitGetTaxonomyFromTree

app = typer.Typer()



class AllTasks(lg.Task):
    """The final node in the Luigi pipeline, creates the DAG."""

    def requires(self):
        output = [
            # V2FastTreeMarkerSplitAnalyseDecoratedFastTree(remove_pct=50)

            # V2FastTreeMarkerSplitCreateMsa(remove_pct=30),

            # BlastGenomeAgainstR207Db(),
            # BlastGetTaxStringAgainstHits(),

            # V2FastTreeMarkerSplitCreateMsa(remove_pct=2),
            # CreateCircosPlot(
            #     ref_gid='GCA_900759445.1',
            #     query_gids=[
            #         'GCF_000735435.1', # planticola
            #         'GCF_001598295.1', # ornitholytica
            #
            #         'GCF_002075345.1',  # citrobacter braaki
            #         'GCF_011064845.1', # freundii
            #                 ],
            #     target_pct=50,
            #     congruence=50
            # ),

        # g__Helicobacter
        # g__Collinsella
        # g__Prevotella
        # V2FastAniInterSpeciesGetValuesForTaxon(taxon='g__Yersinia')
        # CreateCircosPlot(
        #         ref_gid='GCA_014385135.2', # d__Bacteria;p__Coprothermobacterota;c__Coprothermobacteria;o__Coprothermobacterales;f__Coprothermobacteraceae;g__Coprothermobacter;s__Coprothermobacter proteolyticus
        #         query_gids=[
        #             'GCF_014290175.1',  # expected
        #             'GCF_900291955.1',  # o__Oscillospirales;f__Acutalibacteraceae;g__Pseudoruminococcus_A;s__Pseudoruminococcus_A massiliensis
        #         ],
        #         target_pct=50,
        #         congruence=50
        #     ),

            # CreateCircosPlot(
            #     ref_gid='GCF_001812365.1',  # s__Anaerococcus vaginalis
            #     query_gids=[
            #         'GCF_000163295.1',  # s__Anaerococcus vaginalis (rep)
            #         'GCF_000311745.1',  # s__Anaerococcus obesiensis
            #     ],
            #     target_pct=50,
            #     congruence=50
            # ),

            # CreateCircosPlot(
            #     ref_gid='GCA_002937455.1',
            #     query_gids=[
            #         'GCA_002433405.1',
            #         'GCA_002718135.1',
            #     ],
            #     target_pct=50,
            #     congruence=50
            # ),






            FastAniRandomAnalyseResults(target_pct=50)

            # V2FinalGenerateMasterTsv()
            # TODO: 0
            # FastTreeRandomlyRemoveMarkersAnalyseDecorated(target_pct=50, replicate_id=0),
            # V2FastAniRepToClosestRepCreateJobs()

            # V2FastAniInterSpeciesGetValuesForFailedToGenus()

            # FastAniContigSplitZRunForSpecificCase(target_pct=50, query_gid='GCA_900759445.1')

            # V2PplacerMarkerSplitRunPplacer(remove_pct=50),
            # V2PplacerMarkerSplitGetTaxonomyFromTree(remove_pct=50, n_trees_output=8)


        ]

        # for i in range(10):
        #     output.append(FastTreeRandomlyRemoveMarkersAnalyseDecorated(target_pct=50, replicate_id=i))
        # for i in [1, 2, 3, 4, 5, 10, 15, 20, 30, 40]:
        # 1 = howard
        # 2 = bruce
        # 3 = cook
        # 4 = curtin
        # 5 = deakin
        # 10 = hawke
        # 15 = hughes
        # 20 = lyons
        # 30 = mrca002
        # 40 = page

        # for i in [40]:
        #     output.append(V2FastTreeMarkerSplitCreateMsa(remove_pct=i))
        #     output.append(V2PplacerMarkerSplitCreateMsa(remove_pct=i))


        return output

    def complete(self):
        return False


@app.command()
def run():
    typer.echo('Running workflow...')
    lg.build([AllTasks()], workers=1, local_scheduler=True, log_level="WARNING")
    typer.echo('Done.')


@app.command()
def clear():
    if DIR_CACHE is not None:
        delete = typer.confirm(f'Delete directory? {DIR_CACHE}')
        if delete:
            if os.path.isdir(DIR_CACHE):
                shutil.rmtree(DIR_CACHE)
            typer.echo('Done.')
        else:
            typer.echo('Nothing to do.')
    else:
        typer.echo('Cache is not set, nothing to do.')


@app.command()
def rq(queue_name: str):
    typer.echo(f'Running RQ worker: {queue_name}')
    # os.nice(10)

    # Start a worker with a custom name
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        q = Queue(queue_name, connection=conn, default_timeout='60d')
        print(f'Queue size: {len(q)}')
        host_name = socket.gethostname().replace('.ace.uq.edu.au', '')
        worker_name = f'{host_name}-{time.time()}'
        worker = Worker([q], connection=conn, job_monitoring_interval=90, name=worker_name)
        worker.work(burst=True)

    typer.echo('Done.')


@app.command()
def rq_worker_stop(worker_name: str):
    typer.echo(f'Stopping RQ worker: {worker_name}')
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        worker = [x for x in Worker.all(conn) if x.name == worker_name]
        if len(worker) == 0:
            typer.echo('No worker found.')
            return
        worker = worker[0]
        send_shutdown_command(conn, worker)
    typer.echo('Done.')


def main():
    app()


if __name__ == '__main__':
    main()
