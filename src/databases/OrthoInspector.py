import pandas as pd

from databases.Database import Database

class OrthoInspector(Database):
    """OrthoInspector ortholog table between C. elegans and H. sapiens.

    http://lbgi.fr/orthoinspector/dbquery_qfo/data/CSV/62.csv.gz

    Args:
        build_uniprot_list: Whether to build the list of UniProt IDs to search
    """
    def __init__(self, build_uniprot_list=False):
        super().__init__(name="OrthoInspector",
                         filename="orthoinspector",
                         build_uniprot_list=build_uniprot_list)

    def _read_raw(self, build_uniprot_list=False):
        """Returns an ortholog table for OrthoInspector

        Both IDs are provided as Uniprot IDs.

        Returns:
            DataFrame containing the raw orthologs from OrthoInspector
        """
        df = pd.read_csv('data/orthoinspector/orthologs.csv', sep=',',
                         usecols=[0, 1], names=['CE_UNIPROT', 'HS_UNIPROT']) \
                .drop_duplicates()

        # if 'build_uniprot_list' in kwargs and kwargs['build_uniprot_list']:
        if build_uniprot_list:
            self._make_uniprot_list(df)

        return df

    def _perform_worm_mapping(self):
        return pd.merge(self.df, self._get_orthoinspector_uniprot_wb_map(),
                        how='left', on='CE_UNIPROT')

    def _perform_human_mapping(self):
        return pd.merge(self.df, self._get_orthoinspector_uniprot_ensembl_map(),
                        how='left', on='HS_UNIPROT')

    @staticmethod
    def _make_uniprot_list(df):
        """Export list of UniProt IDs to match externally.

        Saves lists of UniProt IDs from worms and humans under

        * `data/orthoinspector/uniprot_list_ce.csv`
        * `data/orthoinspector/uniprot_list_hs.csv`

        respectively.
        """
        uniprot_list_ce = df['CE_UNIPROT'] \
                            .drop_duplicates() \
                            .sort_values()
        uniprot_list_ce. \
            to_csv('data/orthoinspector/uniprot_list_ce.csv', index=False)

        uniprot_list_hs = df['HS_UNIPROT'] \
                            .drop_duplicates() \
                            .sort_values()
        uniprot_list_hs \
            .to_csv('data/orthoinspector/uniprot_list_hs.csv', index=False)

    @staticmethod
    def _get_orthoinspector_uniprot_wb_map():
        """Returns the UniProt to WB ID map for OrthoInspector.

        The first portion is obtained by using the ID mapping tool from UniProt
        available at <http://www.uniprot.org/uploadlists/>. The ones that could
        not be found using this tool are further found by scraping the ID
        history pages.

        Returns:
            A DataFrame with mapping between UniProt and WB IDs
        """
        uniprot_wb_df_1 = pd.read_csv( \
                            'data/orthoinspector/uniprot_wb_map.tsv',
                            sep='\t', header=0, usecols=[0, 1],
                            names=['CE_UNIPROT', 'CE_WB_OLD'])

        # From scraping the history pages
        uniprot_wb_df_2 = pd.read_csv( \
                            'data/orthoinspector/uniprot_wb_map_scraped.csv',
                            sep=',', usecols=[0, 1],
                            names=['CE_UNIPROT', 'CE_WB_OLD'])

        # Combine the two data frames
        uniprot_wb_df = pd.concat([uniprot_wb_df_1, uniprot_wb_df_2], axis=0)

        return uniprot_wb_df

    @staticmethod
    def _get_orthoinspector_uniprot_ensembl_map():
        """Returns the UniProt to Ensembl map for OrthoInspector.

        The first portion is obtained by using the ID mapping tool from UniProt
        available at <http://www.uniprot.org/uploadlists/>. The ones that could
        not be found using this tool are then searched through BioMart tool for
        Ensembl 74 (December 2013), which is closest release after
        OrthoInspector's release (April 2013). The remaining unfound IDs are
        then mapped by scraping the history pages as with
        _get_inparanoid_uniprot_wb_map().

        Returns:
            A DataFrame with mapping between UniProt and Ensembl IDs
        """
        uniprot_ensg_df = pd.read_csv( \
                            'data/orthoinspector/uniprot_ensembl_map.tsv',
                            sep='\t', header=0, usecols=[0, 1],
                            names=['HS_UNIPROT', 'HS_ENSG'])

        # Results putting the "not found" list to Ensembl 74 BioMart:
        #   http://dec2013.archive.ensembl.org/biomart/martview/
        biomart_df_1 = pd.read_csv( \
                        'data/orthoinspector/uniprot_ensembl_map_swissprot.tsv',
                        sep='\t', header=0,
                        names=['HS_ENSG', 'HS_UNIPROT'])
        biomart_df_2 = pd.read_csv( \
                        'data/orthoinspector/uniprot_ensembl_map_trembl.tsv',
                        sep='\t', header=0,
                        names=['HS_ENSG', 'HS_UNIPROT'])
        biomart_df = pd.concat([biomart_df_1, biomart_df_2], axis=0)\
                        .drop_duplicates()

        # Get the scraped map
        scraped_df = pd.read_csv( \
                        'data/orthoinspector/uniprot_ensembl_map_scraped.csv',
                        sep=',', names=['HS_UNIPROT', 'HS_ENSG'])

        # Combine the maps
        uniprot_ensg_df = pd.concat([uniprot_ensg_df, biomart_df, scraped_df],
                                    axis=0) \
                            .drop_duplicates() \
                            .reset_index(drop=True)

        return uniprot_ensg_df
