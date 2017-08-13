import csv
import os.path
import pandas as pd

from databases.Database import Database
from helper.misc import generate_combinations

class InParanoid(Database):
    """InParanoid ortholog table between C. elegans and H. sapiens.

    InParanoid 8.0 (2013-12)
    http://inparanoid.sbc.su.se/download/8.0_current/Orthologs_other_formats/C.elegans/InParanoid.C.elegans-H.sapiens.tgz

    Args:
        build_uniprot_list: Whether to build the list of UniProt IDs to search
    """

    ortholog_file = 'data/inparanoid/orthologs.tsv'

    def __init__(self, build_uniprot_list=False):
        if not os.path.isfile(self.ortholog_file):
            self._make_inparanoid_table(self.ortholog_file)

        super().__init__(name="InParanoid",
                         filename="inparanoid",
                         build_uniprot_list=build_uniprot_list)

    def _read_raw(self, build_uniprot_list=False):
        """Returns an ortholog table for InParanoid

        Both IDs are provided as Uniprot IDs.

        Returns:
            DataFrame containing the raw orthologs from InParanoid
        """
        df = pd.read_csv(self.ortholog_file,
                         sep='\t', header=0, usecols=[0, 1],
                         names=['CE_UNIPROT', 'HS_UNIPROT']) \
                .drop_duplicates()

        if build_uniprot_list:
            self._make_uniprot_list(df)

        return df

    def _perform_worm_mapping(self):
        return pd.merge(self.df, self._get_inparanoid_uniprot_wb_map(),
                        how='left', on='CE_UNIPROT')

    def _perform_human_mapping(self):
        return pd.merge(self.df, self._get_inparanoid_uniprot_ensembl_map(),
                        how='left', on='HS_UNIPROT')

    @staticmethod
    def _make_inparanoid_table(ortholog_file):
        """Creates an ortholog table for InParanoid from the raw SQL table.

        Extracts the C. elegans and H. sapiens orthologs from the raw SQL table.
        Because the orthologs are provided as groupings, combinations are
        generated using generate_combinations(), which uses itertools.product().
        As it is, it current writes to a file first, but it for the future, it
        could be worth directly returning a DataFrame here.

        Args:
            ortholog_file: String of location where the ortholog table will be
            written
        """
        ortholog_list = []
        with open('data/inparanoid/sqltable.C.elegans-H.sapiens') as file:
            reader = csv.reader(file, delimiter='\t')

            current_group = {'group_id': None, 'cele': [], 'hsap': []}

            for row in reader:
                group_id = int(row[0])
                species = row[2]
                uniprot_id = row[4]

                gene_info = uniprot_id

                if group_id != current_group['group_id']:
                    if (len(current_group['cele']) > 0) and \
                        (len(current_group['hsap']) > 0):
                        ortholog_list += generate_combinations(current_group)

                    current_group['group_id'] = group_id
                    current_group['cele'] = []
                    current_group['hsap'] = []

                if species == "C.elegans":
                    current_group['cele'].append(gene_info)
                elif species == "H.sapiens":
                    current_group['hsap'].append(gene_info)

            if (len(current_group['cele']) > 0) \
                and (len(current_group['hsap']) > 0):
                ortholog_list += generate_combinations(current_group)

        with open(ortholog_file, 'w') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(['uniprot_id_cele', 'uniprot_id_hsap'])
            for row in ortholog_list:
                writer.writerow(row)

    @staticmethod
    def _make_uniprot_list(df):
        """Export list of UniProt IDs to match externally.

        Saves lists of UniProt IDs from worms and humans under

        * `data/orthoinspector/uniprot_list_ce.csv`
        * `data/orthoinspector/uniprot_list_hs.csv`

        respectively.
        """
        # Generate lists of UniProt genes to
        uniprot_list_ce = df['CE_UNIPROT'].drop_duplicates().sort_values()
        uniprot_list_ce.to_csv('data/inparanoid/uniprot_list_ce.csv',
                               index=False)

        uniprot_list_hs = df['HS_UNIPROT'].drop_duplicates().sort_values()
        uniprot_list_hs.to_csv('data/inparanoid/uniprot_list_hs.csv',
                               index=False)

    @staticmethod
    def _get_inparanoid_uniprot_wb_map():
        """Returns the UniProt to WB ID map for InParanoid.

        The first portion is obtained by using the ID mapping tool from UniProt
        available at <http://www.uniprot.org/uploadlists/>. The ones that could
        not be found using this tool are further found by scraping the ID
        history pages.

        Returns:
            A DataFrame containing mapping between UniProt and WB IDs.
        """
        uniprot_wb_map_1 = pd.read_csv( \
                            'data/inparanoid/uniprot_wb_map.tsv',
                            sep='\t', header=0, usecols=[0, 1],
                            names=['CE_UNIPROT', 'CE_WB_OLD'])

        # From scraping the history pages
        uniprot_wb_map_2 = pd.read_csv( \
                            'data/inparanoid/uniprot_wb_map_scraped.csv',
                            sep=',', usecols=[0, 1],
                            names=['CE_UNIPROT', 'CE_WB_OLD'])

        # Combine the two data frames
        uniprot_wb_map = pd.concat([uniprot_wb_map_1, uniprot_wb_map_2], axis=0)

        return uniprot_wb_map

    @staticmethod
    def _get_inparanoid_uniprot_ensembl_map():
        """Returns the UniProt to Ensembl map for InParanoid.

        The first portion is obtained by using the ID mapping tool from UniProt
        available at <http://www.uniprot.org/uploadlists/>. The ones that could
        not be found using this tool are then searched through BioMart tool for
        Ensembl 74 (December 2013), which is the closest to InParanoid 8.0's
        release (also December 2013). The remaining unfound IDs are then
        mapped by scraping the history pages as with
        _get_inparanoid_uniprot_wb_map().

        Returns:
            A DataFrame containing mapping between UniProt and Ensembl IDs.
        """
        uniprot_ensg_df = pd.read_csv( \
                            'data/inparanoid/uniprot_ensembl_map.tsv',
                            sep='\t', header=0, usecols=[0, 1],
                            names=['HS_UNIPROT', 'HS_ENSG'])

        # Results putting the "not found" list to Ensembl 74 BioMart:
        #   http://dec2013.archive.ensembl.org/biomart/martview/
        biomart_df_1 = pd.read_csv( \
                        'data/inparanoid/uniprot_ensembl_map_swissprot.tsv',
                        sep='\t', header=0, names=['HS_ENSG', 'HS_UNIPROT'])
        biomart_df_2 = pd.read_csv( \
                        'data/inparanoid/uniprot_ensembl_map_trembl.tsv',
                        sep='\t', header=0, names=['HS_ENSG', 'HS_UNIPROT'])
        biomart_df = pd.concat([biomart_df_1, biomart_df_2], axis=0) \
                        .drop_duplicates()

        # Get the scraped map
        scraped_df = pd.read_csv( \
                        'data/inparanoid/uniprot_ensembl_map_scraped.csv',
                        sep=',', names=['HS_UNIPROT', 'HS_ENSG'])

        # Combine the maps
        uniprot_ensg_df = pd.concat([uniprot_ensg_df, biomart_df, scraped_df],
                                    axis=0) \
                            .drop_duplicates().reset_index(drop=True)

        return uniprot_ensg_df
