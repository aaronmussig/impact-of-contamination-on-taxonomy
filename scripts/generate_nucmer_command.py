from collections import defaultdict

from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.util.paths import get_gid_root
import os
import pandas as pd


GIDS = ['GCA_016712675.1']


R95_SP_TO_R207_REP = {
's__Ruminococcus_H bromii_A': 's__UBA1417 sp004555625',
's__GCA-900066135 sp900543575': 's__Limivivens sp900543575',
's__UBA4871 sp900554535': 's__Caproicibacterium sp900554535',
's__Dorea sp900553355': 's__Merdimonas sp900553355',
's__UMGS1670 sp900553995': 's__Metalachnospira sp900553995',
's__Agathobacter sp900550845': 's__Agathobacter sp905209075',
's__Dialister sp900555245': 's__Dialister invisus',
's__Agathobacter sp900557055': 's__Agathobacter faecis',
's__Dorea sp900550865': 's__Dorea_A sp900550865',
's__Dorea sp900066765': 's__Lachnoclostridium_B sp900066765',
's__Dorea longicatena_B': 's__Dorea_A longicatena_B',
's__Faecalimonas nexilis': 's__Faecalimonas phoceensis',
's__Raoultella planticola': 's__Klebsiella planticola',
's__Raoultella ornithinolytica': 's__Raoultella ornithinolytica',
's__Raoultella electrica': 's__Klebsiella electrica',
's__Klebsiella_A grimontii': 's__Klebsiella grimontii',
's__Serratia_B rubidaea_A': 's__Serratia_B rhizosphaerae',
's__159R sp004346745': 's__Acerihabitans sp004346745',
's__Raoultella sp003752615': 's__Klebsiella sp003752615',
's__Raoultella sp004342285': 's__Klebsiella sp004342285',
's__Klebsiella_A michiganensis_B': 's__Klebsiella spallanzanii',
's__Klebsiella_A oxytoca_B': 's__Klebsiella huaxiensis',
's__Raoultella terrigena': 's__Klebsiella terrigena',
's__Raoultella sp002270295': 's__Klebsiella sp002270295',
's__Klebsiella_A michiganensis': 's__Klebsiella michiganensis',
's__Citrobacter sp004684345': 's__Citrobacter tructae',
's__Serratia plymuthica_A': 's__Serratia inhibens',
's__Pantoea sp002389975': 's__Pantoea piersonii',
's__Acidianus_A hospitalis': 's__Acidianus hospitalis',
's__Acidianus_A manzaensis': 's__Acidianus manzaensis',
's__QEFN01 sp003086555': 's__Stygiolobus sp003086555',
's__Sulfolobus_B sp001316085': 's__Sulfuracidifex tepidarius',
's__Acidianus_A sulfidivorans': 's__Acidianus sulfidivorans',
's__Thermofilum pendens': 's__Thermofilum_B pendens_A',
's__Sulfolobus_B metallicus': 's__Sulfuracidifex metallicus',
    's__Dorea longicatena': 's__Dorea_A longicatena',
    's__Mogibacterium sp002299625': 's__VUNA01 sp002299625',
    's__Muricomes sp900604355': 's__Luxibacter massiliensis',
    's__UBA2856 sp900555005': 's__Porcincola sp900555005',
}
def main():

    df_meta = GtdbMetadataR207().output().read_cached()
    df_reps = df_meta[df_meta['gtdb_representative'] == 't']

    for gid in GIDS:

        gid_root = get_gid_root(gid)
        gid_fna = os.path.join(gid_root, f'{gid}.fna')

        cmd = [
            'mummer2circos',
            '-l',
            '-r',
            gid_fna,
            '-q',
        ]

        path_assign = os.path.join(gid_root, 'gunc_r95', 'gunc_output', f'{gid}.contig_assignments.tsv')
        df_contig = pd.read_csv(path_assign, sep='\t')
        df_contig = df_contig[df_contig['tax_level'] == 'species']

        d_assignment_counts = defaultdict(list)
        for _, row in df_contig.iterrows():
            d_assignment_counts[row['contig']].append(
                (row['assignment'], row['count_of_genes_assigned'])
            )

        taxa_to_keep = set()
        d_contig_to_counts = dict()
        for contig, assignments in d_assignment_counts.items():
            # Find the proportion of this contig
            total = sum(count for _, count in assignments)
            cur_counts = dict()
            for taxon, count in assignments:
                taxon_pct = count / total * 100
                cur_counts[taxon] = taxon_pct
                # if taxon_pct >= 50:
                #     taxa_to_keep.add(taxon)
            d_contig_to_counts[contig] = cur_counts

        taxa_to_keep = set()
        for _, row in df_contig.sort_values(by='count_of_genes_assigned', ascending=False).iterrows():
            taxa_to_keep.add(row['assignment'])

        for _, row in df_contig.iterrows():
            species = row['assignment']
            adj_species = R95_SP_TO_R207_REP.get(species, species)
            # print(adj_species)
            try:
                cur_rep = df_reps[df_reps['species'] == adj_species].index[0]
                cur_rep_path = os.path.join(get_gid_root(cur_rep), f'{cur_rep}.fna')
                cmd.append(cur_rep_path)
                print(cur_rep, adj_species)
            except Exception:
                print('?', adj_species)


        print(gid)
        print(' '.join(cmd))
        print()

        break

    pass




if __name__ == '__main__':
    main()