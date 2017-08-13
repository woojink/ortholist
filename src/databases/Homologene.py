import csv
import os.path
import pandas as pd
import re

from databases.Database import Database
from helper.misc import generate_combinations, tidy_split

class Homologene(Database):
    """Homologene ortholog table between C. elegans and H. sapiens.

    Build 68, 2014-05-06
    ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data

    Args:
        build_entrez_list: Whether to build the list of Entrez IDs to search
    """

    ortholog_file = 'data/homologene/orthologs.tsv'

    def __init__(self, build_entrez_list=False):
        if not os.path.isfile(self.ortholog_file):
            self._make_homologene_table(self.ortholog_file)

        super().__init__(name="Homologene",
                         filename="homologene",
                         build_entrez_list=build_entrez_list)

    def _read_raw(self, build_entrez_list=False):
        """Returns an ortholog table for Homologene

        Both IDs are provided as Entrez IDs.

        Returns:
            DataFrame containing the raw orthologs from Homologene
        """
        df = pd.read_csv(self.ortholog_file, sep='\t', header=0, usecols=[0, 1],
                             names=['CE_ENTREZ', 'HS_ENTREZ']) \
                .drop_duplicates()

        if build_entrez_list:
            self._make_entrez_list(df)

        return df

    def _perform_worm_mapping(self):
        return pd.merge(self.df, self._get_homologene_entrez_wb_map(),
                          how='left', on='CE_ENTREZ')

    def _perform_human_mapping(self):
        return pd.merge(self.df, self._get_homologene_entrez_ensembl_map(),
                          how='left', on='HS_ENTREZ')

    @staticmethod
    def _make_homologene_table(ortholog_file):
        """Creates an ortholog table for Homologene from the raw table.

        Extracts the C. elegans and H. sapiens orthologs from the raw table.
        Because the orthologs are provided as groupings, combinations are
        generated using generate_combinations(), which uses itertools.product().
        As it is, it current writes to a file first, but it for the future, it
        would be a good idea directly returning a DataFrame here.

        Args:
            ortholog_file (str): Location where the ortholog table will be written
        """
        ortholog_list = []
        with open('data/homologene/homologene.tsv') as file:
            reader = csv.reader(file, delimiter='\t')

            current_group = {'group_id': None, 'cele': [], 'hsap': []}

            for row in reader:
                group_id = int(row[0])
                taxonomy_id = int(row[1])
                entrez_id = int(row[2])

                gene_info = entrez_id

                if group_id != current_group['group_id']:
                    if (len(current_group['cele']) > 0) and \
                            (len(current_group['hsap']) > 0):
                        ortholog_list += generate_combinations(current_group)

                    current_group['group_id'] = group_id
                    current_group['cele'] = []
                    current_group['hsap'] = []

                if taxonomy_id == 6239:
                    current_group['cele'].append(gene_info)
                elif taxonomy_id == 9606:
                    current_group['hsap'].append(gene_info)

            if (len(current_group['cele']) > 0) \
                    and (len(current_group['hsap']) > 0):
                ortholog_list += generate_combinations(current_group)


        with open(ortholog_file, 'w') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(['entrez_id_cele', 'entrez_id_hsap'])
            for row in ortholog_list:
                writer.writerow(row)

    @staticmethod
    def _make_entrez_list(df):
        """Export list of UniProt IDs to match externally.

        Saves lists of UniProt IDs from worms and humans under

        * `data/homologene/entrez_list_ce.csv`
        * `data/homologene/entrez_list_hs.csv`

        respectively.
        """
        entrez_list_ce = df['CE_ENTREZ']\
                            .drop_duplicates() \
                            .sort_values()
        entrez_list_ce.to_csv('data/homologene/entrez_list_ce.csv',
                                index=False)

        entrez_list_hs = df['HS_ENTREZ'] \
                            .drop_duplicates() \
                            .sort_values()
        entrez_list_hs.to_csv('data/homologene/entrez_list_hs.csv',
                                index=False)

    @staticmethod
    def _get_homologene_entrez_wb_map():
        """Returns the Entrez to WB ID map for Homologene.

        The first portion is obtained by using gene info table from
        <ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz>,
        downloaded on 2017-01-28. The ones that could not be found here are
        further found by scraping the ID history pages.

        Returns:
            A DataFrame containing mapping between Entrez and WB IDs
        """
        entrez_wb_df = pd.read_csv( \
                        'data/entrez/Caenorhabditis_elegans.gene_info.gz',
                        sep='\t', header=0, usecols=[1, 5],
                        names=['CE_ENTREZ', 'CE_WB_OLD'])

        # Pick out WB ID entries and separate with commas
        entrez_wb_df['CE_WB_OLD'] = entrez_wb_df['CE_WB_OLD'] \
                                    .apply(Homologene._get_wb_id_for_entrez)

        # Split comma-separated values in one row into multiple rows
        entrez_wb_df = tidy_split(entrez_wb_df, 'CE_WB_OLD', sep=',')

        # Load the rest from the scraped results
        scraped_df = pd.read_csv( \
                        'data/homologene/entrez_wb_map_scraped.csv', 
                        sep=',',
                        names=['CE_ENTREZ', 'CE_WB_OLD'])

        entrez_wb_df = pd.concat([entrez_wb_df, scraped_df], axis=0)

        return entrez_wb_df

    @staticmethod
    def _get_homologene_entrez_ensembl_map():
        """Returns the Entrez to Ensembl map for Homologene.

        The first portion is obtained by using gene info table from
        <ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz>,
        downloaded on 2017-01-28. The ones that could not be found here are
        further found by scraping the ID history pages.

        Returns:
            A DataFrame containing mapping between Entrez and Ensembl IDs
        """
        entrez_ensg_df = pd.read_csv( \
                            'data/entrez/Homo_sapiens.gene_info.gz',
                            sep='\t', header=0, usecols=[1, 5],
                            names=['HS_ENTREZ', 'HS_ENSG'])

        # Pick out WB ID entries and separate with commas
        entrez_ensg_df['HS_ENSG'] = entrez_ensg_df['HS_ENSG'] \
                                    .apply(Homologene._get_ensg_id_for_entrez)

        # Split comma-separated values in one row into multiple rows
        entrez_ensg_df = tidy_split(entrez_ensg_df, 'HS_ENSG', sep=',')

        # Load the rest from the scraped results
        scraped_df = pd.read_csv( \
                        'data/homologene/entrez_ensembl_map_scraped.csv',
                        sep=',', names=['HS_ENTREZ', 'HS_ENSG'])

        entrez_ensg_df = pd.concat([entrez_ensg_df, scraped_df], axis=0)

        return entrez_ensg_df

    @staticmethod
    def _get_wb_id_for_entrez(entry):
        """Finds all the WB ID entries and returns them separated with commas."""
        id_set = re.findall(r'WormBase:(WBGene[0-9]{8})', entry)
        if len(id_set):
            return ','.join(id_set)
        return None

    @staticmethod
    def _get_ensg_id_for_entrez(entry):
        """Finds all the ENSG entries and returns them separated with commas."""
        id_set = re.findall(r'ENSG[0-9]{11}', entry)
        if len(id_set):
            return ','.join(id_set)
        return None
