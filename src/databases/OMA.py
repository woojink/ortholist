import csv
import gzip
import pandas as pd

from databases.Database import Database

class OMA(Database):
    """OMA ortholog table between C. elegans and H. sapiens.
    """
    def __init__(self):
        super().__init__(name="OMA", filename="oma")

    def _read_raw(self):
        """Returns an ortholog table for OMA

        Human IDs are provided as ENSG, but worm IDs are provided as Wormpep.

        Returns:
            A DataFrame containing the raw orthologs from OMA
        """
        return pd.read_csv('data/oma/orthologs.tsv', sep='\t', header=None,
                usecols=[0, 1], names=['CE_WORMPEP', 'HS_ENSG'])

    def _perform_worm_mapping(self):
        return pd.merge(self.df, self._get_oma_wb_map(), left_on='CE_WORMPEP',
                right_index=True).drop('CE_WORMPEP', axis=1)

    @staticmethod
    def _get_oma_wb_map():
        """Returns OMA to WB ID mapping

        Wormpep data from Wormbase:
        ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/protein/c_elegans.WS235.wormpep_package.tar.gz

        OMA to Wormpep data from OMA:
        http://omabrowser.org/All/oma-wormbase.txt.gz

        Returns:
            A DataFrame with mapping between OMA ID to Wormbase ID
        """
        wp_to_wb = {}
        with gzip.open('data/wormbase/wormpep.table235.gz', 'rt') as file:
            reader = csv.reader(file, delimiter="\t")
            for row in reader:
                wp_to_wb[row[1]] = row[2]

        oma_to_wb = {}
        with gzip.open('data/oma/oma-wormbase.txt.gz', 'rt') as file:
            reader = csv.reader(file, delimiter="\t")
            next(reader)
            next(reader)
            for oma_id, wp_id in reader:
                if not wp_id.startswith('CE'):
                    continue
                oma_to_wb[oma_id] = wp_to_wb[wp_id]

        oma_to_wb = pd.DataFrame.from_dict(oma_to_wb, orient='index')
        oma_to_wb.columns = ['CE_WB_OLD']

        return oma_to_wb
