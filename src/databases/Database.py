import csv
import gzip
import pandas as pd

from collections import defaultdict

from helper.misc import generate_combinations
from helper.wb_map import get_ce_wb_updated

class Database(object):
    """An ortholog database.

    Args:
        name: Name of the database
    """

    def __init__(self, name, filename, **kwargs):
        self.name = name
        self.filename = filename
        self.df = self._read_raw(**kwargs)
        self.df = self._perform_worm_mapping()
        self.df = self._perform_human_mapping()
        self.df = self._process_wb_changes(self.df)
        self.df = self.df.drop_duplicates()

    def _process_wb_changes(self, df):
        """Processes the database with WormBase ID updates.

        Returns:
            A DataFrame with WormBase IDs mapped to either current IDs or None if
            deprecated. Also includes the old IDs and comments if changed.
        """
        # Deal with WB ID changes
        df = pd.concat([df, get_ce_wb_updated(df)], axis=1) \
                    .sort_values(['CE_WB_CURRENT', 'HS_ENSG'])

        # Return the rearranged database
        return df[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD', 'CE_WB_COMMENT']]

    def get_df(self):
        """Returns the final database

        The columns are as follows:
            * CE_WB_CURRENT: Current WormBase ID of the worm gene
            * HS_ENSG: Ensembl ID of the human gene
            * CE_WB_OLD: Referenced WormBase ID from the database
            * CE_WB_COMMENT: Reason (if any) of the WormBase ID change

        Returns:
            The final processed database DataFrame
        """
        return self.df

    def _read_raw(self):
        """Reads the raw database.

        Raises:
            NotImplementedError: Abstract method here, should be implemented for
            by the implementing class.
        """
        raise NotImplementedError()

    def _perform_worm_mapping(self):
        return self.df

    def _perform_human_mapping(self):
        return self.df
