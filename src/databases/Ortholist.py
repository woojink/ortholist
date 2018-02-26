import pandas as pd

from databases.Database import Database

class Ortholist(Database):
    """Legacy Ortholist 1.0 database"""
    def __init__(self):
        super().__init__(name="Ortholist 1.0", filename="ortholist_1")

    def _process_wb_changes(self):
        # Override because these are deprecated by default
        return self.df

    def _read_raw(self):
        """Returns the raw Ortholist database

        Returns:
            Raw DataFrame of Ensembl Compara release 87-89
        """
        source = 'data/ortholist/legacy.csv'

        # Read the ortholog list
        df = pd.read_csv(source, header=0, usecols=[0, 1]).drop_duplicates()

        return df
