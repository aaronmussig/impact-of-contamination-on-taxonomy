from typing import Literal

from workflow.config import R207_AR53_HMM, R207_BAC120_HMM


def get_genome_domain_from_markers(results_markers) -> Literal['d__Bacteria', 'd__Archaea']:
    """Determine domain of User genomes based on identified marker genes."""

    single_copy_hits = {**results_markers['muq'], **results_markers['unq']}

    arc_count = len(set(single_copy_hits).intersection(set(R207_AR53_HMM)))
    bac_count = len(set(single_copy_hits).intersection(set(R207_BAC120_HMM)))

    arc_aa_pct = arc_count * 100 / len(R207_AR53_HMM)
    bac_aa_pct = bac_count * 100 / len(R207_BAC120_HMM)

    if bac_aa_pct >= arc_aa_pct:
        return 'd__Bacteria'
    else:
        return 'd__Archaea'
