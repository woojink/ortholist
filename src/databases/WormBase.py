import csv
import gzip
import pandas as pd

from collections import defaultdict

from helper.wb_map import get_ce_wb_updated

class WormBase(object):
    """Table for WormBase ID to common name mapping.

        WormBase gene ID and InterPro domain tables (2016-09-15) from:
        <ftp://ftp.wormbase.org/pub/wormbase/releases/WS255/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz>
        <ftp://ftp.wormbase.org/pub/wormbase/releases/WS255/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS255.protein_domains.tsv>

        This table contains WormBase IDs, common names, as well as locus ID
        information. Ahringer RNAi clone locations (WS239) are read, updated,
        and left-joined with the WormBase table to provide non-ortholog
        C. elegans information.
    """

    def __init__(self):
        self.name = "WormBase"
        self.filename = "wormbase"
        self.db = self._get_wormbase_table()

    def get_df(self):
        """Returns the consolidated WormBase table

        Returns:
            A DataFrame containing mapping information between WB IDs,
            common names, locus IDs, and Ahringer RNAi clone locations.
        """
        return self.db

    @staticmethod
    def _get_wormbase_table():
        """Returns the consolidated WormBase table

        Returns:
            A DataFrame containing mapping information between WB IDs,
            common names, locus IDs, and Ahringer RNAi clone locations.
        """
        wb_df = pd.read_csv( \
                    "data/wormbase/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz",
                    sep=',', compression='gzip', header=None, usecols=[1, 2, 3],
                    names=['CE_WB_CURRENT', 'COMMON_NAME', 'LOCUS_ID'])

        # Read Ahringer locations, mapped to WS239
        ahringer_df = pd.read_csv("data/ahringer/locations_ws239.csv", sep=',',
                                  header=None, usecols=[0,1],
                                  names=["AHRINGER_LOC", 'CE_WB_OLD']) \
                        .groupby('CE_WB_OLD')['AHRINGER_LOC'] \
                        .apply(lambda x: '|'.join(x)) \
                        .reset_index()

        # Deal with WB ID changes
        ahringer_df = pd.concat([ahringer_df, get_ce_wb_updated(ahringer_df)], \
                            axis=1) \
                        .sort_values(['CE_WB_CURRENT']).reset_index(drop=True)

        # Read InterPro domains
        ip_dict = defaultdict(set)
        with gzip.open('data/wormbase/c_elegans.PRJNA13758.WS255.protein_domains.tsv.gz', 'rt') as file:
            reader = csv.reader(file, delimiter="\t")
            for line in reader:
                wb_id = line[0]
                domains = line[3:]

                # InterPro domains are separated by '|' for each WormBase ID
                for domain in domains:
                    ip_dict[wb_id].add(domain)
        interpro_df = pd.DataFrame( \
                        [(k, '|'.join(sorted(v))) for k,v in ip_dict.items()],
                        columns=['CE_WB_CURRENT', 'INTERPRO_DOM'])

        # Left join using the WormBase gene ID table
        wb_df = pd.merge(wb_df, ahringer_df, how='left', on='CE_WB_CURRENT')
        wb_df = pd.merge(wb_df, interpro_df, how='left', on='CE_WB_CURRENT')

        wb_df = wb_df[['CE_WB_CURRENT', 'COMMON_NAME', 'LOCUS_ID', 
                        'AHRINGER_LOC', 'INTERPRO_DOM']]

        return wb_df
