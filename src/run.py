#!/usr/bin/env python3
"""Script to process OrthoList 2 and write to various files
"""
import csv
import gzip

from collections import defaultdict

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database

from databases.EnsemblCompara import EnsemblCompara
from databases.Homologene import Homologene
from databases.InParanoid import InParanoid
from databases.OMA import OMA
from databases.OrthoInspector import OrthoInspector
from databases.Ortholist import Ortholist
from databases.OrthoMCL import OrthoMCL
from databases.WormBase import WormBase


ENSEMBL_LOCATION = 'data/ensembl/89/ensembl_annotations.tsv.gz'
OMIM_LOCATION = 'data/omim/OMIMDATA_2018-02-07.csv'


def write_to_csv(df, filename, gzip=False):
    """Write the DataFrame to a CSV

    Args:
        df: DataFrame to write
        filename: String for the filename to write to
        gzip: Boolean of whether to compress with gzip or not
    """
    if not gzip:
        df.to_csv('results/{filename}.csv'.format(filename=filename),
                  index=False)
    else:
        df.to_csv('results/{filename}.csv.gz'.format(filename=filename),
                  index=False, compression='gzip')

def get_ensembl_annotations():
    """Retrieve SMART, GO, and HGNC information from Ensembl 89

    Returns:
        A DataFrame containing SMART, GO, and HGNC annotations from Ensembl 89
    """
    ensembl_data = defaultdict(lambda: defaultdict(set))
    with gzip.open(ENSEMBL_LOCATION, 'rt') as file:
        reader = csv.reader(file, delimiter="\t")
        for ensg, smart, go_terms, _, hgnc in reader:
            if smart:
                ensembl_data[ensg]['SMART'].add(smart)
            if go_terms:
                ensembl_data[ensg]['GO'].add(go_terms)
            if 'HGNC' not in ensembl_data[ensg] and hgnc:
                ensembl_data[ensg]['HGNC'] = hgnc
    ensembl_df = pd.DataFrame.from_dict(ensembl_data, orient='index')
    ensembl_df.reset_index(inplace=True)
    ensembl_df.rename(columns={'index': 'HS_ENSG'}, inplace=True)

    return ensembl_df


def get_omim_annotations():
    """Retrieve OMIM annotations

    Returns:
        A DataFrame containing OMIM annotations
    """
    omim_data = defaultdict(lambda: defaultdict(list))
    with open(OMIM_LOCATION, 'r') as file:
        reader = csv.reader(file)
        next(reader)

        for ensg, gene, phenotype in reader:
            omim_data[ensg]['OMIM_GENES'] = set(gene.split())
            omim_data[ensg]['OMIM_PHENOTYPES'] = \
                set([x.strip() for x in phenotype.split('|')])
    omim_df = pd.DataFrame.from_dict(omim_data, orient='index')
    omim_df.reset_index(inplace=True)
    omim_df.rename(columns={'index': 'HS_ENSG'}, inplace=True)

    return omim_df

if __name__ == "__main__":
    ####################
    # Process the orthologs
    ####################
    print('\nProcessing the orthologs...')
    COMPARA = EnsemblCompara()
    HOMOLOGENE = Homologene()
    INPARANOID = InParanoid()
    OMA_DF = OMA()
    ORTHOINSPECTOR = OrthoInspector()
    ORTHOMCL = OrthoMCL()
    WORMBASE = WormBase()

    ORTHOLOG_DATABASES = [COMPARA, HOMOLOGENE, INPARANOID,
                          OMA_DF, ORTHOINSPECTOR, ORTHOMCL]
    ALL_DATABASES = ORTHOLOG_DATABASES + [WORMBASE]

    ####################
    # Write to CSV
    ####################
    print('\nWriting to CSV...')
    for db in ALL_DATABASES:
        df, name, filename = db.get_df(), db.name, db.filename
        print("    Writing {name}".format(name=db.name))
        write_to_csv(df, filename)
    print("Done!")

    ####################
    # Write to Excel
    ####################
    print('\nWriting to Excel...')
    WRITER = pd.ExcelWriter('results/results.xlsx')
    for db in ALL_DATABASES:
        df, name, filename = db.get_df(), db.name, db.filename
        print("    Preparing {name}".format(name=db.name))
        df.to_excel(WRITER, name, index=False)
    WRITER.save()
    print("Done!")

    ####################
    # Make combined database
    ####################
    print('\nProcessing combined database...')
    COMBINED_DF = pd.DataFrame()
    for db in ORTHOLOG_DATABASES:
        COMBINED_DF = COMBINED_DF.append(db.get_df()).drop_duplicates()
    print("    Writing combined CSV")
    write_to_csv(COMBINED_DF, "combined")

    # Export unique IDs
    print("    Writing unique IDs")
    write_to_csv(COMBINED_DF['HS_ENSG'] \
        .drop_duplicates().sort_values(), "unique_ensg")
    write_to_csv(COMBINED_DF['CE_WB_CURRENT'] \
        .drop_duplicates().sort_values(), "unique_wb")
    print("Done")

    ####################
    # Master table
    ####################
    print('\nPreparing master table...')

    ALL_PAIRS = defaultdict(set)
    for db in ORTHOLOG_DATABASES:
        df, name = db.get_df(), db.name
        nn_lines = df['CE_WB_CURRENT'].notnull() & df['HS_ENSG'].notnull()
        for ce, hs in df[['CE_WB_CURRENT', 'HS_ENSG']][nn_lines].values:
            ALL_PAIRS[(ce, hs)].add(name)

    ## Create a consolidated pair list with list of databases and a score
    ##  MASTER_TUPLES contains the following:
    ##   (worm gene, human gene, database list, score [number of databases])
    MASTER_TUPLES = [(k[0], k[1], sorted(list(v)), len(v)) \
                        for k, v in ALL_PAIRS.items()]
    MASTER_DF = pd.DataFrame(MASTER_TUPLES,
                             columns=['CE_WB_CURRENT', 'HS_ENSG',
                                      'Databases', 'Score'])

    ## List of overlap present in Ensembl Compara 89
    ENSEMBL89_ENSG = pd.read_csv('data/ensembl/89/ensg_list.csv',
                                 header=0, names=["HS_ENSG"])

    ## Throw away ENSG IDs not present in Ensembl Compara 89 per Dan
    MASTER_DF = pd.merge(MASTER_DF, ENSEMBL89_ENSG, how="inner", on="HS_ENSG")

    ## Fetch legacy orthologs and append
    ORTHOLIST_DF = Ortholist().get_df()
    ORTHOLIST_DF['Databases'] = [['Legacy Ortholist'] for _ in range(len(ORTHOLIST_DF))]
    ORTHOLIST_DF['Score'] = 0
    MASTER_DF = MASTER_DF.append(ORTHOLIST_DF)

    ## Add information from WormBase db (common name, Ahringer location, etc.)
    MASTER_DF = pd.merge(MASTER_DF, WORMBASE.get_df(),
                         how='left', on='CE_WB_CURRENT')

    ## Add information from Ensembl 89 annotations (SMART, GO, HGNC name)
    MASTER_DF = pd.merge(MASTER_DF, get_ensembl_annotations(),
                         how='left', on='HS_ENSG')

    ## Add information from OMIM annotations
    MASTER_DF = pd.merge(MASTER_DF, get_omim_annotations(),
                         how='left', on='HS_ENSG')

    ## Join lists into pipe-separated strings
    MASTER_DF['Databases'] = MASTER_DF['Databases'] \
        .apply(lambda x: '|'.join(sorted(list(x))) \
                if isinstance(x, list) else None)
    MASTER_DF['SMART'] = MASTER_DF['SMART'] \
        .apply(lambda x: '|'.join(sorted(list(x))) \
                if isinstance(x, set) else None)
    MASTER_DF['GO'] = MASTER_DF['GO'] \
        .apply(lambda x: '|'.join(sorted(list(x), key=lambda x: x.lower())) \
                if isinstance(x, set) else None)
    MASTER_DF['OMIM_GENES'] = MASTER_DF['OMIM_GENES'] \
        .apply(lambda x: '|'.join(sorted(list(x))) \
                if isinstance(x, set) else None)
    MASTER_DF['OMIM_PHENOTYPES'] = MASTER_DF['OMIM_PHENOTYPES'] \
           .apply(lambda x: '|'.join(sorted(list(x), key=lambda x: x.lower())) \
                if isinstance(x, set) else None)

    ## Write to CSV
    print("    Writing to CSV")
    write_to_csv(MASTER_DF, "master", gzip=True)

    ## Write to database
    print('\nConnecting to database...')
    ENGINE = create_engine('mysql://root:@localhost/ortholist')
    if not database_exists(ENGINE.url):
        create_database(ENGINE.url)
    MASTER_DF.to_sql(
        name='ortholist',
        con=ENGINE,
        if_exists='replace',
    )

    ## Write to Excel, merge rows with multi-indexing
    print("    Writing to Excel")
    WRITER = pd.ExcelWriter('results/master.xlsx')
    MASTER_DF.set_index(['CE_WB_CURRENT']) \
             .to_excel(WRITER)
    WRITER.save()
    print('Done!')
