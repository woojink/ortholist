import pandas as pd

from databases.Database import Database

class EnsemblCompara(Database):
    """Ensembl Compara database release 87 through 89

    Ensembl 87 BioMart (2016-12)
    http://dec2016.archive.ensembl.org/biomart/martview

    Ensembl 88 BioMart (2017-03)
    http://mar2017.archive.ensembl.org/biomart/martview

    Ensembl 89 BioMart (2017-05)
    http://may2017.archive.ensembl.org/biomart/martview
    """
    _VERSIONS = [87, 88, 89]

    def __init__(self):
        super().__init__(name="Ensembl Compara 87-89", filename="compara87-89")

    def _read_raw(self):
        """Returns the raw Ensembl Compara database

        Returns:
            Raw DataFrame of Ensembl Compara release 87-89
        """
        source = 'data/ensembl/{version}/orthologs.tsv'
        df = pd.DataFrame()

        for version in self._VERSIONS:
            # Read the ortholog list
            df= df.append(pd.read_csv(source.format(version=version),
                            sep='\t', header=0, usecols=[0, 2],
                            names=['CE_WB_OLD', 'HS_ENSG'])).drop_duplicates()
        return df
