#!/usr/bin/env python3

import csv
import gzip
import os.path
import re
import subprocess

import pandas as pd

from collections import defaultdict
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database

from databases.EnsemblCompara import EnsemblCompara
from databases.Homologene import Homologene
from databases.InParanoid import InParanoid
from databases.OMA import OMA
from databases.OrthoInspector import OrthoInspector
from databases.OrthoMCL import OrthoMCL
from databases.WormBase import WormBase

from helper.misc import generate_combinations, tidy_split
from helper.wb_map import get_ce_wb_updated

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
    with gzip.open('data/ensembl/89/ensembl_annotations.tsv.gz', 'rt') as file:
        reader = csv.reader(file, delimiter="\t")
        for ensg, smart, go, _, hgnc in reader:
            if smart:
                ensembl_data[ensg]['SMART'].add(smart)
            if go:
                ensembl_data[ensg]['GO'].add(go)
            if 'HGNC' not in ensembl_data[ensg] and hgnc:
                ensembl_data[ensg]['HGNC'] = hgnc
    ensembl_df = pd.DataFrame.from_dict(ensembl_data, orient='index')
    ensembl_df.reset_index(inplace=True)
    ensembl_df.rename(columns={'index': 'HS_ENSG'}, inplace=True)

    return ensembl_df

if __name__ == "__main__":
    ####################
    # Process the orthologs
    ####################
    print('\nProcessing the orthologs...')
    COMPARA = EnsemblCompara()
    HOMOLOGENE = Homologene()
    INPARANOID = InParanoid()
    OMA = OMA()
    ORTHOINSPECTOR = OrthoInspector()
    ORTHOMCL = OrthoMCL()
    WORMBASE = WormBase()
    
    ORTHOLOG_DATABASES = [COMPARA, HOMOLOGENE, INPARANOID,
                          OMA, ORTHOINSPECTOR, ORTHOMCL]
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
    writer = pd.ExcelWriter('results/results.xlsx')
    for db in ALL_DATABASES:
        df, name, filename = db.get_df(), db.name, db.filename
        print("    Preparing {name}".format(name=db.name))
        df.to_excel(writer, name, index=False)
    writer.save()
    print("Done!")

    ####################
    # Make combined database
    ####################
    print('\nProcessing combined database...')
    combined_df = pd.DataFrame()
    for db in ORTHOLOG_DATABASES:
        combined_df = combined_df.append(db.get_df()).drop_duplicates()
    print("    Writing combined CSV")
    write_to_csv(combined_df, "combined")

    # Export unique IDs
    print("    Writing unique IDs")
    write_to_csv(combined_df['HS_ENSG'].drop_duplicates(), "unique_ensg")
    write_to_csv(combined_df['CE_WB_CURRENT'].drop_duplicates(), "unique_wb")
    print("Done")

    ####################
    # Master table
    ####################
    print('\nPreparing master table...')

    all_pairs = defaultdict(set)
    for db in ORTHOLOG_DATABASES:
        df, name = db.get_df(), db.name
        nn_lines = df['CE_WB_CURRENT'].notnull() & df['HS_ENSG'].notnull()
        for ce, hs in df[['CE_WB_CURRENT','HS_ENSG']][nn_lines].values:
            all_pairs[(ce, hs)].add(name)

    ## Create a consolidated pair list with list of databases and a score
    ##  master_tuples contains the following:
    ##   (worm gene, human gene, database list, score [number of databases])
    master_tuples = [(k[0], k[1], sorted(list(v)), len(v)) \
                        for k,v in all_pairs.items()]
    master_df = pd.DataFrame(master_tuples,
                    columns=['CE_WB_CURRENT', 'HS_ENSG', 'Databases', 'Score'])

    ## List of overlap present in Ensembl Compara 89
    ENSEMBL89_ENSG = pd.read_csv('data/ensembl/89/ensg_list.csv',
                                    header=0, names=["HS_ENSG"])

    ## Throw away ENSG IDs not present in Ensembl Compara 89 per Dan
    master_df = pd.merge(master_df, ENSEMBL89_ENSG, how="inner", on="HS_ENSG")

    ## Add information from WormBase db (common name, Ahringer location, etc.)
    master_df = pd.merge(master_df, WORMBASE.get_df(),
                            how='left', on='CE_WB_CURRENT')

    ## Add information from Ensembl 89 annotations (SMART, GO, HGNC name)
    master_df = pd.merge(master_df, get_ensembl_annotations(),
                            how='left', on='HS_ENSG')

    ## Join lists into pipe-separated strings
    master_df['Databases'] = master_df['Databases'] \
        .apply(lambda x: '|'.join(sorted(list(x))) \
                if isinstance(x, list) else None)
    master_df['SMART'] = master_df['SMART'] \
        .apply(lambda x: '|'.join(sorted(list(x))) \
                if isinstance(x, set) else None)
    master_df['GO'] = master_df['GO'] \
        .apply(lambda x: '|'.join(sorted(list(x), key=lambda x: x.lower())) \
                if isinstance(x, set) else None)

    ## Write to CSV
    print("    Writing to CSV")
    write_to_csv(master_df, "master", gzip=True)

    ## Write to Excel, merge rows with multi-indexing
    print("    Writing to Excel")
    writer = pd.ExcelWriter('results/master.xlsx')
    master_df.set_index(['CE_WB_CURRENT', 'COMMON_NAME', 'LOCUS_ID',
                            'AHRINGER_LOC', 'INTERPRO_DOM', 'HS_ENSG']) \
             .to_excel(writer)
    writer.save()
    print('Done!')
