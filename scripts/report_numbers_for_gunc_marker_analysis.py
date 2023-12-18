from tqdm import tqdm

from workflow.fasttree_marker_analysis.gunc_d_analyse_decorated import FastTreeMarkerAnalyseDecorated
from workflow.fasttree_marker_analysis.random_d_analyse_decorated import FastTreeMarkerAnalyseDecoratedRandom
import pandas as pd

def parse(obj):
    batch_id = obj.batch_id
    path = obj.output().path

    df = pd.read_csv(path, sep='\t')

    all_gids = set(df['gid'])
    gids_to_check = set(df[df['tax_result'] == 'check']['gid'])
    gids_congruent = set(df[df['tax_result'] == 'congruent']['gid'])
    gids_same_msa = set(df[df['classification'] == 'same_msa']['gid'])
    gids_changed_domain = set(df[df['classification'] == 'changed_domain']['gid'])
    gids_no_markers = set(df[df['classification'] == 'no_markers']['gid'])

    gids_run_on = set(df[df['classification'] == 'run_on']['gid'])
    assert(len(gids_to_check.intersection(gids_congruent).intersection(gids_same_msa).intersection(gids_changed_domain).intersection(gids_no_markers)) == 0)
    assert(len(gids_to_check.intersection(gids_run_on)) == len(gids_to_check))
    assert(len(gids_to_check.union(gids_congruent).union(gids_same_msa).union(gids_changed_domain).union(gids_no_markers).union(gids_run_on)) == len(all_gids))

    gids_congruent = gids_run_on - gids_to_check

    print(f'Batch {batch_id}')
    print('\n'.join(sorted(gids_to_check)))


    return {
        'batch': batch_id,
        'no_markers': len(gids_no_markers),
        'same_msa': len(gids_same_msa),
        'changed_domain': len(gids_changed_domain),
        'congruent': len(gids_congruent),
        'incongruent': len(gids_to_check),
    }


def get_gids(obj: FastTreeMarkerAnalyseDecorated):
    path = obj.output().path

    df = pd.read_csv(path, sep='\t')

    all_fail = set()

    all_gids = set(df['gid'])
    gids_to_check = set(df[df['tax_result'] == 'check']['gid'])
    gids_congruent = set(df[df['tax_result'] == 'congruent']['gid'])
    gids_same_msa = set(df[df['classification'] == 'same_msa']['gid'])
    gids_changed_domain = set(df[df['classification'] == 'changed_domain']['gid'])
    gids_no_markers = set(df[df['classification'] == 'no_markers']['gid'])

    gids_run_on = set(df[df['classification'] == 'run_on']['gid'])
    assert (len(gids_to_check.intersection(gids_congruent).intersection(gids_same_msa).intersection(
        gids_changed_domain).intersection(gids_no_markers)) == 0)
    assert (len(gids_to_check.intersection(gids_run_on)) == len(gids_to_check))
    assert (len(gids_to_check.union(gids_congruent).union(gids_same_msa).union(gids_changed_domain).union(
        gids_no_markers).union(gids_run_on)) == len(all_gids))

    gids_congruent = gids_run_on - gids_to_check

    print(F'All gids: {len(all_gids)}')
    print(F'No markers: {len(gids_no_markers)}')
    print(F'Same MSA: {len(gids_same_msa)}')
    print(F'Changed domain: {len(gids_changed_domain)}')
    print(F'Congruent: {len(gids_congruent)}')
    print(F'Incongruent: {len(gids_to_check)}')
    # print(gids_to_check)
    print()

    return gids_to_check


def main():
    rows = list()
    for i in tqdm(range(0, 10)):  # 20
        rows.append(parse(FastTreeMarkerAnalyseDecoratedRandom(target_pct=50, batch_id=i)))

    df = pd.DataFrame(rows)
    #
    # for col in ['no_markers', 'same_msa', 'changed_domain', 'congruent', 'incongruent']:
    #     print(f'{col} = {df[col].mean():,}')
    #
    # df.to_csv('/tmp/stats2.tsv', sep='\t', index=False)

    rows = set()
    rows.update(get_gids(FastTreeMarkerAnalyseDecorated(target_pct=50, congruence=50)))
    # rows.update(get_gids(FastTreeMarkerAnalyseDecorated(target_pct=50, congruence=50)))

    print('\n'.join(sorted(rows)))
    print(len(rows))

    return


if __name__ == '__main__':
    main()