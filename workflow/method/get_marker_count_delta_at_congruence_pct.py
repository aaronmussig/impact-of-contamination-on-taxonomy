from tqdm import tqdm

from workflow.config import R207_BAC120_HMM
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.genome import Genome
from workflow.model.gunc_model import GuncRefDb


def _get_base_markers_for_gid(gid):
    genome = Genome(gid)
    d_marker_hits = genome.get_marker_hits()
    out = set()
    out.update(set(d_marker_hits['unq'].keys()))
    out.update(set(d_marker_hits['muq'].keys()))
    return frozenset(out)


def run_on_gid(gid, congruence, pct, max_css, source):

    domain_markers = frozenset(R207_BAC120_HMM)

    base_markers = _get_base_markers_for_gid(gid)
    base_markers = base_markers.intersection(domain_markers)

    genome = Genome(gid)
    contigs, pct_removed = genome.get_gunc_contigs_where_removed_equals_x_pct_genome_removed(
        pct=pct,
        max_congruence=congruence,
        max_css=max_css,
        source=source
    )
    contigs = frozenset(contigs)


    new_markers = genome.get_marker_hits(contigs)
    new_unq_hits = {**new_markers['muq'], **new_markers['unq']}
    new_unq_hits = set(new_unq_hits.keys())
    new_unq_hits = new_unq_hits.intersection(domain_markers)

    markers_lost = base_markers - new_unq_hits
    markers_gained = new_unq_hits - base_markers

    return base_markers, markers_lost, markers_gained, pct_removed



def main():
    gid_set = {'GCA_012964215.1', 'GCA_018822835.1'}

    import pandas as pd

    df_css = AggregateMaxCssLevelMerged().output().maybe_read_cached()
    rows = list()
    for gid in tqdm(gid_set):

        row = df_css.loc[gid]

        if row['source'] == 'progenomes':
            source = GuncRefDb.PRO
        else:
            source = GuncRefDb.GTDB

        base_markers, markers_lost, markers_gained, pct_removed = run_on_gid(gid, congruence=50, pct=50, max_css=row['taxonomic_level'], source=source)

        rows.append({
            'gid': gid,
            'base_markers': len(base_markers),
            'markers_lost': len(markers_lost),
            'markers_gained': len(markers_gained),
            'pct_removed': pct_removed
        })

    df = pd.DataFrame(rows)
    df.sort_values(by='gid', inplace=True)

    return

if __name__ == '__main__':


    main()

