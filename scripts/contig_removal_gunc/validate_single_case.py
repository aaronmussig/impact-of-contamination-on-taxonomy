from gunc.get_scores import read_diamond_output
from gunc.gunc import run_gunc
from magna.util.disk import copy_file

from workflow.config import PATH_GUNC_GTDB_REF_DB
from workflow.gunc_helper.aggregate_diamond_output import GuncGtdbR95AggregateDiamondOutput
from workflow.util.paths import get_gid_root_cutoff
import os
import pandas as pd


def generate_data():
    gid = 'GCA_002495045.1'
    cutoffs = [0, 1, 84, 86, 91]

    target_dir = '/srv/home/uqamussi/tmp/gunc_validate_single_case'
    target_genome_dir = os.path.join(target_dir, 'genomes')
    target_output_dir = os.path.join(target_dir, 'output')
    os.makedirs(target_genome_dir, exist_ok=True)
    os.makedirs(target_output_dir, exist_ok=True)

    for cutoff in cutoffs:
        gid_root = get_gid_root_cutoff(gid, cutoff)
        faa_path = os.path.join(gid_root, 'prodigal', f'{gid}.faa')

        new_path = os.path.join(target_genome_dir, f'{gid}_{cutoff}.faa')
        copy_file(faa_path, new_path, checksum=True)

    cmd = ['gunc', 'run', '--db_file', PATH_GUNC_GTDB_REF_DB,
           '--file_suffix', '.faa', '--gene_calls',
           '--input_dir', target_genome_dir, '--threads',
           str(90), '--detailed_output',
           '--contig_taxonomy_output', '--out_dir', target_output_dir,
           '--temp_dir', '/tmp']

    print(' '.join(cmd))


def validate():
    gid = 'GCA_002495045.1'
    cutoffs = [0, 1, 84, 86, 91]
    genes_called = [1646, 1594, 1588, 1581, 1570]
    target_dir = '/srv/home/uqamussi/tmp/gunc_validate_single_case'
    target_genome_dir = os.path.join(target_dir, 'genomes')
    output_dir = os.path.join(target_dir, 'output')
    diamond_dir = os.path.join(output_dir, 'diamond_output')

    path_diamond_df = GuncGtdbR95AggregateDiamondOutput().output().path
    df_diamond: pd.DataFrame = pd.read_hdf(path_diamond_df, where=f'gid == "{gid}"')
    df_diamond.drop(columns='gid', inplace=True)
    df_diamond.to_csv(f'/tmp/{gid}.diamond.gtdb_95.out', sep='\t', index=False, header=False)


    d_kv = dict()
    with open('/srv/home/uqamussi/tmp/gunc_validate_single_case/output/diamond_output/GCA_002495045.1_0.diamond.gtdb_95.out') as f:
        for line in f.readlines():
            cols = line.strip().split('\t')
            if cols[0] not in d_kv:
                d_kv[cols[0]] = line
            else:
                if line != d_kv[cols[0]]:
                    print('???')

    with open(f'/tmp/{gid}.diamond.gtdb_95.out') as f:
        for line in f.readlines():
            cols = line.strip().split('\t')
            if cols[0] not in d_kv:
                d_kv[cols[0]] = line
            else:
                other = d_kv[cols[0]]
                if line[:50] != other[:50]:
                    print('???')


    d_kv = dict()

    for i, cutoff in enumerate(cutoffs):
        diamond_path = os.path.join(diamond_dir, f'{gid}_{cutoff}.diamond.gtdb_95.out')

        df_gunc = read_diamond_output(diamond_path)

        df_max_css = run_gunc(diamond_outfiles=[f'/tmp/{gid}.diamond.gtdb_95.out'],
                              genes_called={f'{gid}': genes_called[i]},
                              out_dir='/tmp',
                              sensitive=False,
                              detailed_output=True,
                              db='gtdb_95',
                              min_mapped_genes=11,
                              use_species_level=False,
                              contig_taxonomy_output=True)



        print()


        # with open(diamond_path, 'r') as f:
        #     for line in f:
        #         cols = line.strip().split('\t')
        #
        #         if line[0] in d_kv:
        #             if d_kv[cols[0]] != line:
        #                 print('???')
        #         else:
        #             d_kv[(cols[0])] = line


def main():

    # generate_data()
    validate()
    pass


if __name__ == '__main__':
    main()