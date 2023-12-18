import dendropy
from gtdblib.util.bio.accession import canonical_gid
from tqdm import tqdm

from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.fasttree_marker_split.e_decorate_fasttree import FastTreeMarkerSplitDecorateFastTree
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged

import os

from workflow.util.log import log


class ArbFilterFile:

    def __init__(self, path: str):
        self.path = path
        self.columns = list()

    def write(self):
        lines = [
            'AUTODETECT\t"BEGIN"',
            '',
            'BEGIN\t"BEGIN*"',
        ]

        for col in self.columns:
            lines.extend([
                '',
                f'MATCH\t"{col}\=*"',
                '\tSRT "*\=="',
            ])
            if col == 'db_name':
                lines.append(f'\tWRITE "name"')
            else:
                lines.append(f'\tWRITE "{col}"')

        lines.extend([
            '',
            f'SEQUENCEAFTER\t"{self.columns[-1]}*"',
            'SEQUENCESRT\t"*\=="',
            'SEQUENCEEND\t"END"',
            '',
            'END\t"END"',
            ''
        ])

        with open(self.path, 'w') as f:
            f.write('\n'.join(lines))


class ArbMetadataFile:

    def __init__(self, path):
        self.path = path

    def write(self):
        return


def create_dummy_tree(path: str, taxa):
    with open(path, 'w') as f:
        taxa = ','.join(taxa)
        f.write(f'({taxa});')
    return


def main():

    out_dir = '/tmp/export'
    os.makedirs(out_dir, exist_ok=True)

    log('Reading tree')
    tree = FastTreeMarkerSplitDecorateFastTree().output().read()
    fail_gids_canonical = {canonical_gid(x.label.replace('TEST_', '').replace('_C', '')) for x in tree.taxon_namespace if x.label.startswith('TEST_')}

    # Update the taxon names
    for taxon in tree.taxon_namespace:
        if taxon.label.startswith('TEST_'):
            if taxon.label.endswith('_C'):
                taxon.label = canonical_gid(taxon.label[5:-2]) + "_C"
            else:
                taxon.label = canonical_gid(taxon.label[5:])
        else:
            taxon.label = canonical_gid(taxon.label)

    tree.write(path=os.path.join(out_dir, 'fasttree_marker_split_halves.tree'), schema='newick', suppress_rooting=True)

    gids_to_write = set()
    for gid in fail_gids_canonical:
        gids_to_write.add(gid)
        gids_to_write.add(gid + '_C')

    log('Creating dummy tree')
    create_dummy_tree(os.path.join(out_dir, 'failed_genome_null.tree'), gids_to_write)

    arb_filter = ArbFilterFile(os.path.join(out_dir, 'arb_filter.ift'))

    log('Reading css')
    df_css = AggregateMaxCssLevelMerged().output().maybe_read_cached()

    log('Reading meta')
    df_meta = GtdbMetadataR207().output().maybe_read_cached()

    log('Merging dfs')
    df_merged = df_meta.merge(df_css, left_index=True, right_index=True)

    # Create the arb filter
    log('Creating arb filter')
    arb_filter.columns.append('db_name')
    for col in df_merged.columns:
        if col == 'pass.GUNC':
            col = 'pass_gunc'
        arb_filter.columns.append(col)
    arb_filter.write()

    # Write out the data for the arb metadata file
    log('Writing arb metadata file')
    with open(os.path.join(out_dir, 'arb_metadata.txt'), 'w') as f:
        for gid, row in tqdm(df_merged.iterrows(), total=len(df_merged)):
            gid = canonical_gid(gid)
            if gid in fail_gids_canonical:
                for suffix in ('', '_C'):
                    f.write('BEGIN\n')
                    f.write(f'db_name={gid}{suffix}\n')
                    for col, val in row.items():
                        if col == 'pass.GUNC':
                            col = 'pass_gunc'
                        f.write(f'{col}={val}\n')
                    f.write('aligned_seq=A\n')
                    f.write('END\n\n')

    return


if __name__ == '__main__':
    main()
