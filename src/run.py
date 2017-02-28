#!/usr/bin/env python3

import csv
import gzip
import os.path
import re
import subprocess

import pandas as pd

from helper.misc import generate_combinations, tidy_split
from helper.wb_map import get_ce_wb_updated

from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database

####################
# OrthoMCL
####################
def get_orthomcl():
    """Returns the OrthoMCL ortholog table between C. elegans and H. sapiens.

    The C. elegans side is provided as ENSP IDs instead of ENSG IDs, so the IDs are mapped using
    Wormbase release 56, which OrthoMCL version 5 uses.

    Version 5 (2015-07-23):
    http://orthomcl.org/common/downloads/release-5/pairs/orthologs.txt.gz

    Returns:
        pandas.DataFrame: Data frame containing the orthologs from OrthoMCL
    """

    # Read the ortholog list
    orthomcl = pd.read_csv('data/orthomcl/orthologs.csv', names=['CE_WB_OLD', 'HS_ENSP'])

    # Convert ENSP to ENSG
    orthomcl = pd.merge(orthomcl, _get_ensembl_56_ensp_ensg_map(),
                        left_on='HS_ENSP', right_index=True)

    # Deal with WB ID changes
    orthomcl = pd.concat([orthomcl, get_ce_wb_updated(orthomcl)], axis=1) \
                 .sort_values(['CE_WB_CURRENT', 'HS_ENSG'])

    # Rearrange the columns
    orthomcl = orthomcl[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD', 'CE_WB_COMMENT']]

    return orthomcl

def _get_ensembl_56_ensp_ensg_map():
    """Return ENSP to ENSG mapping from Ensembl release 56

    Because there convinient direct mapping from ENSP to ENSG, various tables are used to
    map between the different stages.

    Data from
    ftp://ftp.ensembl.org/pub/release-56/mysql/homo_sapiens_core_56_37a/

    The schema is described in
        data/ensembl/56/homo_sapiens_core_56_37a.sql.gz
    """
    ensp_ensg_df = pd.read_csv("data/ensembl/56/translation_stable_id.txt.gz", sep='\t',
                               compression='gzip', header=None, usecols=[0, 1],
                               names=["translation_id", "ENSP"])

    translation_to_transcript = pd.read_csv("data/ensembl/56/translation.txt.gz", sep='\t',
                                            compression='gzip', header=None, usecols=[0, 1],
                                            names=["translation_id", "transcript_id"])
    transcript_to_gene = pd.read_csv("data/ensembl/56/transcript.txt.gz", sep='\t',
                                     compression='gzip', header=None, usecols=[0, 1],
                                     names=["transcript_id", "gene_id"])
    gene_id_to_gene = pd.read_csv("data/ensembl/56/gene_stable_id.txt.gz", sep='\t',
                                  compression='gzip', header=None, usecols=[0, 1],
                                  names=["gene_id", 'HS_ENSG'])

    ensp_ensg_df = ensp_ensg_df.set_index("translation_id") \
        .join(translation_to_transcript.set_index("translation_id")) \
        .set_index("transcript_id").join(transcript_to_gene.set_index("transcript_id")) \
        .set_index("gene_id").join(gene_id_to_gene.set_index("gene_id")).set_index("ENSP")

    return ensp_ensg_df

####################
# OMA
####################
def get_oma():
    """Returns the OMA ortholog table between C. elegans and H. sapiens.

    OMA, May 2016 release
    http://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=CAEEL&p2=HUMAN&p3=EnsemblGene

    Returns:
        pandas.DataFrame: DataFrame containing the orthologs from OMA
    """

    # Read the ortholog list
    oma = pd.read_csv('data/oma/orthologs.tsv', sep='\t', header=None, usecols=[0, 1],
                      names=['CE_WORMPEP', 'HS_ENSG'])

    # Convert OMA IDs to Wormbase
    oma = pd.merge(oma, _get_oma_wb_map(),
                   left_on='CE_WORMPEP', right_index=True).drop('CE_WORMPEP', axis=1)

    # Deal with WB ID changes
    oma = pd.concat([oma, get_ce_wb_updated(oma)], axis=1) \
            .sort_values(['CE_WB_CURRENT', 'HS_ENSG'])

    # Rearrange the columns
    oma = oma[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD', 'CE_WB_COMMENT']]

    return oma

def _get_oma_wb_map():
    """Returns OMA to WB ID mapping

    Wormpep data from Wormbase:
    ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/protein/c_elegans.WS235.wormpep_package.tar.gz

    OMA to Wormpep data from OMA:
    http://omabrowser.org/All/oma-wormbase.txt.gz

    Returns:
        pandas.DataFrame: Mapping data from OMA ID to Wormbase ID
    """
    wp_to_wb = {}
    with open('data/wormbase/wormpep.table235') as file:
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

####################
# Ensembl Compara
####################
def get_compara():
    """Returns the Ensembl Compara ortholog table between C. elegans and H. sapiens.

    Ensembl 87 BioMart (2016-12)
    http://dec2016.archive.ensembl.org/biomart/martview

    Returns:
        pandas.DataFrame: DataFrame containing the orthologs from Ensembl Compara
    """

    # Read the ortholog list
    compara = pd.read_csv('data/ensembl/87/orthologs.tsv', sep='\t', header=0, usecols=[0, 2],
                          names=['CE_WB_OLD', 'HS_ENSG'])

    # Deal with WB ID changes
    compara = pd.concat([compara, get_ce_wb_updated(compara)], axis=1) \
                .sort_values(['CE_WB_CURRENT', 'HS_ENSG'])

    # Rearrange the columns
    compara = compara[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD', 'CE_WB_COMMENT']]

    return compara

####################
# InParanoid
####################
def get_inparanoid(uniprot_map_built=True):
    """Returns the InParanoid ortholog table between C. elegans and H. sapiens.

    InParanoid 8.0 (2013-12)
    http://inparanoid.sbc.su.se/download/8.0_current/Orthologs_other_formats/C.elegans/InParanoid.C.elegans-H.sapiens.tgz

    Returns:
        pandas.DataFrame: DataFrame containing the orthologs from InParanoid
    """
    ortholog_file = 'data/inparanoid/orthologs.tsv'

    if not os.path.isfile(ortholog_file):
        _make_inparanoid_table(ortholog_file)

    # Read the ortholog list
    inparanoid = pd.read_csv(ortholog_file, sep='\t', header=0, usecols=[0, 1],
                             names=['CE_UNIPROT', 'HS_UNIPROT'])

    if not uniprot_map_built:
        # Generate lists of UniProt genes to
        inparanoid_uniprot_list_ce = inparanoid['CE_UNIPROT'].drop_duplicates().sort_values()
        inparanoid_uniprot_list_ce.to_csv('data/inparanoid/uniprot_list_ce.csv', index=False)
        print("Number of unique worm UniProt IDs: {}".format(len(inparanoid_uniprot_list_ce)))

        inparanoid_uniprot_list_hs = inparanoid['HS_UNIPROT'].drop_duplicates().sort_values()
        inparanoid_uniprot_list_hs.to_csv('data/inparanoid/uniprot_list_hs.csv', index=False)
        print("Number of unique worm UniProt IDs: {}".format(len(inparanoid_uniprot_list_hs)))

        # Use the list of UniProt IDs to get WB IDs from
        #   http://www.uniprot.org/uploadlists/

    # Convert C.elegans UniProt to WB ID
    inparanoid = pd.merge(inparanoid, _get_inparanoid_uniprot_wb_map(),
                          how='left', on='CE_UNIPROT')

    # Convert H. sapiens UniProt to Ensembl
    inparanoid = pd.merge(inparanoid, _get_inparanoid_uniprot_ensembl_map(),
                          how='left', on='HS_UNIPROT')

    # Deal with WB ID changes
    inparanoid = pd.concat([inparanoid, get_ce_wb_updated(inparanoid)], axis=1) \
            .sort_values(['CE_WB_CURRENT', 'HS_ENSG'])

    # Rearrange the columns
    inparanoid = inparanoid[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD',
                             'CE_WB_COMMENT', 'CE_UNIPROT', 'HS_UNIPROT']]

    return inparanoid

def _make_inparanoid_table(ortholog_file):
    """Creates an ortholog table for InParanoid from the raw SQL table.

    Extracts the C. elegans and H. sapiens orthologs from the raw SQL table.
    Because the orthologs are provided as groupings, combinations are generated
    using generate_combinations(), which uses itertools.product(). As it is, it
    current writes to a file first, but it for the future, it could be worth
    directly returning a DataFrame here.

    Args:
        ortholog_file (str): Location where the ortholog table will be written
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
                if (len(current_group['cele']) > 0) and (len(current_group['hsap']) > 0):
                    ortholog_list += generate_combinations(current_group)

                current_group['group_id'] = group_id
                current_group['cele'] = []
                current_group['hsap'] = []

            if species == "C.elegans":
                current_group['cele'].append(gene_info)
            elif species == "H.sapiens":
                current_group['hsap'].append(gene_info)

        if (len(current_group['cele']) > 0) and (len(current_group['hsap']) > 0):
            ortholog_list += generate_combinations(current_group)

    with open(ortholog_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['uniprot_id_cele', 'uniprot_id_hsap'])
        for row in ortholog_list:
            writer.writerow(row)

def _get_inparanoid_uniprot_wb_map():
    """Returns the UniProt to WB ID map for InParanoid.

    The first portion is obtained by using the ID mapping tool from UniProt available at
    <http://www.uniprot.org/uploadlists/>. The ones that could not be found using this tool
    are further found by scraping the ID history pages.

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between UniProt and WB IDs.
    """
    uniprot_wb_map_1 = pd.read_csv('data/inparanoid/uniprot_wb_map.tsv', sep='\t', header=0,
                                   usecols=[0, 1], names=['CE_UNIPROT', 'CE_WB_OLD'])

    # From scraping the history pages
    uniprot_wb_map_2 = pd.read_csv('data/inparanoid/uniprot_wb_map_scraped.tsv', sep=',',
                                   usecols=[0, 1], names=['CE_UNIPROT', 'CE_WB_OLD'])

    # Combine the two data frames
    uniprot_wb_map = pd.concat([uniprot_wb_map_1, uniprot_wb_map_2], axis=0)

    return uniprot_wb_map

def _get_inparanoid_uniprot_ensembl_map(uniprot_map_built=True):
    """Returns the UniProt to Ensembl map for InParanoid.

    The first portion is obtained by using the ID mapping tool from UniProt available at
    <http://www.uniprot.org/uploadlists/>. The ones that could not be found using this tool
    are then searched through BioMart tool for Ensembl 74 (December 2013), which is the closest
    to InParanoid 8.0's release (also December 2013). The remaining unfound IDs are then
    mapped by scraping the history pages as with _get_inparanoid_uniprot_wb_map().

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between
            UniProt and Ensembl IDs.
    """
    uniprot_ensg_df = pd.read_csv('data/inparanoid/uniprot_ensembl_map.tsv', sep='\t',
                                  header=0, usecols=[0, 1], names=['HS_UNIPROT', 'HS_ENSG'])

    # Results putting the "not found" list to Ensembl 74 BioMart:
    #   http://dec2013.archive.ensembl.org/biomart/martview/
    biomart_df_1 = pd.read_csv('data/inparanoid/uniprot_ensembl_map_swissprot.tsv',
                               sep='\t', header=0, names=['HS_ENSG', 'HS_UNIPROT'])
    biomart_df_2 = pd.read_csv('data/inparanoid/uniprot_ensembl_map_trembl.tsv',
                               sep='\t', header=0, names=['HS_ENSG', 'HS_UNIPROT'])
    biomart_df = pd.concat([biomart_df_1, biomart_df_2], axis=0).drop_duplicates()

    if not uniprot_map_built:
        # Compare the "not found" list from UniProt with the BioMart set generated above
        # (135 are still missing), export the difference
        not_found_set = set(pd.read_csv('data/inparanoid/uniprot_ensembl_map_not_found.tsv',
                                        sep='\t', header=0, names=['HS_UNIPROT'])['HS_UNIPROT'])
        biomart_set = set(biomart_df['HS_UNIPROT'])

        pd.Series(list(not_found_set.difference(biomart_set))).sort_values()\
            .to_csv('data/inparanoid/uniprot_ensembl_map_to_scrape.csv', index=False)

    # Get the scraped map
    scraped_df = pd.read_csv('data/inparanoid/uniprot_ensembl_map_scraped.csv',
                             sep=',', names=['HS_UNIPROT', 'HS_ENSG'])

    # Combine the maps
    uniprot_ensg_df = pd.concat([uniprot_ensg_df, biomart_df, scraped_df], axis=0) \
                        .drop_duplicates().reset_index(drop=True)

    return uniprot_ensg_df

####################
# OrthoInspector
####################
def get_orthoinspector(preprocessed=True, uniprot_map_built=True):
    """OrthoInspector

    http://lbgi.fr/orthoinspector/dbquery_qfo/data/CSV/62.csv.gz

    Returns:
        pandas.DataFrame: DataFrame containing the orthologs from OrthoInspector
    """
    if not preprocessed:
        print('Downloading and preprocessing OrthoInspector orthologs...\n')
        subprocess.Popen(['./preprocess.sh'], cwd='data/orthoinspector', shell=True)

    # Read the orthologs
    orthoinspector = pd.read_csv('data/orthoinspector/orthologs.csv', sep=',',
                                 usecols=[0, 1], names=['CE_UNIPROT', 'HS_UNIPROT'])

    if not uniprot_map_built:
        orthoinspector_uniprot_list_ce = orthoinspector['CE_UNIPROT'].drop_duplicates().sort_values()
        orthoinspector_uniprot_list_ce. \
            to_csv('data/orthoinspector/uniprot_list_ce.csv', index=False)
        print("Number of unique worm UniProt IDs: {}".format(len(orthoinspector_uniprot_list_ce)))

        orthoinspector_uniprot_list_hs = orthoinspector['HS_UNIPROT'].drop_duplicates().sort_values()
        orthoinspector_uniprot_list_hs \
            .to_csv('data/orthoinspector/uniprot_list_hs.csv', index=False)
        print("Number of unique human UniProt IDs: {}".format(len(orthoinspector_uniprot_list_hs)))

    # Convert C.elegans UniProt to WB ID
    orthoinspector = pd.merge(orthoinspector, _get_orthoinspector_uniprot_wb_map(),
                              how='left', on='CE_UNIPROT')

    # Convert H. sapiens Uniprot to WB ID
    orthoinspector = pd.merge(orthoinspector,
                              _get_orthoinspector_uniprot_ensembl_map(uniprot_map_built),
                              how='left', on='HS_UNIPROT')

    # Deal with WB ID changes
    orthoinspector = pd.concat([orthoinspector, get_ce_wb_updated(orthoinspector)], axis=1) \
            .sort_values(['CE_WB_CURRENT', 'HS_ENSG']).reset_index(drop=True)

    # Rearrange the columns
    orthoinspector = orthoinspector[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD',
                                     'CE_WB_COMMENT', 'CE_UNIPROT', 'HS_UNIPROT']]

    return orthoinspector

def _get_orthoinspector_uniprot_wb_map():
    """Returns the UniProt to WB ID map for OrthoInspector.

    The first portion is obtained by using the ID mapping tool from UniProt available at
    <http://www.uniprot.org/uploadlists/>. The ones that could not be found using this tool
    are further found by scraping the ID history pages.

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between UniProt and WB IDs.
    """
    uniprot_wb_df_1 = pd.read_csv('data/orthoinspector/uniprot_wb_map.tsv', sep='\t', header=0,
                                  usecols=[0, 1], names=['CE_UNIPROT', 'CE_WB_OLD'])

    # From scraping the history pages
    uniprot_wb_df_2 = pd.read_csv('data/orthoinspector/uniprot_wb_map_scraped.csv', sep=',',
                                  usecols=[0, 1], names=['CE_UNIPROT', 'CE_WB_OLD'])

    # Combine the two data frames
    uniprot_wb_df = pd.concat([uniprot_wb_df_1, uniprot_wb_df_2], axis=0)

    return uniprot_wb_df

def _get_orthoinspector_uniprot_ensembl_map(uniprot_map_built=True):
    """Returns the UniProt to Ensembl map for OrthoInspector.

    The first portion is obtained by using the ID mapping tool from UniProt available at
    <http://www.uniprot.org/uploadlists/>. The ones that could not be found using this tool
    are then searched through BioMart tool for Ensembl 74 (December 2013), which is closest
    release after OrthoInspector's release (April 2013). The remaining unfound IDs are then
    mapped by scraping the history pages as with _get_inparanoid_uniprot_wb_map().

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between
            UniProt and Ensembl IDs.
    """
    uniprot_ensg_df = pd.read_csv('data/orthoinspector/uniprot_ensembl_map.tsv', sep='\t',
                                  header=0, usecols=[0, 1], names=['HS_UNIPROT', 'HS_ENSG'])

    # Results putting the "not found" list to Ensembl 74 BioMart:
    #   http://dec2013.archive.ensembl.org/biomart/martview/
    biomart_df_1 = pd.read_csv('data/orthoinspector/uniprot_ensembl_map_swissprot.tsv',
                               sep='\t', header=0, names=['HS_ENSG', 'HS_UNIPROT'])
    biomart_df_2 = pd.read_csv('data/orthoinspector/uniprot_ensembl_map_trembl.tsv',
                               sep='\t', header=0, names=['HS_ENSG', 'HS_UNIPROT'])
    biomart_df = pd.concat([biomart_df_1, biomart_df_2], axis=0).drop_duplicates()

    if not uniprot_map_built:
        # Compare the "not found" list from UniProt with the BioMart set generated above
        # (13 are still missing), export the difference
        not_found_set = set(pd.read_csv('data/orthoinspector/uniprot_ensembl_map_not_found.tsv',
                                        sep='\t', header=0, names=['HS_UNIPROT'])['HS_UNIPROT'])
        biomart_set = set(biomart_df['HS_UNIPROT'])

        pd.Series(list(not_found_set.difference(biomart_set))).sort_values() \
            .to_csv('data/orthoinspector/uniprot_ensembl_map_to_scrape.csv', index=False)

    # Get the scraped map
    scraped_df = pd.read_csv('data/orthoinspector/uniprot_ensembl_map_scraped.csv',
                             sep=',', names=['HS_UNIPROT', 'HS_ENSG'])

    # Combine the maps
    uniprot_ensg_df = pd.concat([uniprot_ensg_df, biomart_df, scraped_df], axis=0) \
                        .drop_duplicates().reset_index(drop=True)

    return uniprot_ensg_df

####################
# Homologene
####################
def get_homologene(preprocessed=True):
    """Homologene

    Build 68, 2014-05-06
    ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data

    Returns:
        pandas.DataFrame: DataFrame containing the orthologs from Homologene
    """
    ortholog_file = 'data/homologene/orthologs.tsv'

    if not preprocessed:
        print('Preprocessing Homologene orthologs...')
        subprocess.Popen(['./preprocess.sh'], cwd='data/homologene', shell=True)

    if not os.path.isfile(ortholog_file):
        _make_homologene_table(ortholog_file)

    # Read the ortholog list
    homologene = pd.read_csv(ortholog_file, sep='\t', header=0, usecols=[0, 1],
                             names=['CE_ENTREZ', 'HS_ENTREZ'])

    homologene_entrez_list_ce = homologene['CE_ENTREZ'].drop_duplicates().sort_values()
    homologene_entrez_list_ce.to_csv('data/homologene/entrez_list_ce.csv', index=False)
    print("Number of unique worm Entrez IDs: {}".format(len(homologene_entrez_list_ce)))

    homologene_entrez_list_hs = homologene['HS_ENTREZ'].drop_duplicates().sort_values()
    homologene_entrez_list_hs.to_csv('data/homologene/entrez_list_hs.csv', index=False)
    print("Number of unique human Entrez IDs: {}".format(len(homologene_entrez_list_hs)))

    # Convert C.elegans Entrez to WB ID
    homologene = pd.merge(homologene, _get_homologene_entrez_wb_map(homologene_entrez_list_ce),
                          how='left', on='CE_ENTREZ')

    # Convert H. Sapiens Entrez to WB ID
    homologene = pd.merge(homologene, _get_homologene_entrez_ensembl_map(homologene_entrez_list_hs),
                          how='left', on='HS_ENTREZ')

    # Deal with WB ID changes
    homologene = pd.concat([homologene, get_ce_wb_updated(homologene)], axis=1) \
            .sort_values(['CE_WB_CURRENT', 'HS_ENSG']).reset_index(drop=True)

    # Rearrange the columns
    homologene = homologene[['CE_WB_CURRENT', 'HS_ENSG', 'CE_WB_OLD',
                             'CE_WB_COMMENT', 'CE_ENTREZ', 'HS_ENTREZ']]

    return homologene

def _make_homologene_table(ortholog_file):
    """Creates an ortholog table for Homologene from the raw table.

    Extracts the C. elegans and H. sapiens orthologs from the raw table.
    Because the orthologs are provided as groupings, combinations are generated
    using generate_combinations(), which uses itertools.product(). As it is, it
    current writes to a file first, but it for the future, it would be a good idea
    directly returning a DataFrame here.

    Args:
        ortholog_file (str): Location where the ortholog table will be written
    """
    ortholog_list = []
    with open('data/homologene/homologene.txt') as file:
        reader = csv.reader(file, delimiter='\t')

        current_group = {'group_id': None, 'cele': [], 'hsap': []}

        for row in reader:
            group_id = int(row[0])
            taxonomy_id = int(row[1])
            entrez_id = int(row[2])

            gene_info = entrez_id

            if group_id != current_group['group_id']:
                if (len(current_group['cele']) > 0) and (len(current_group['hsap']) > 0):
                    ortholog_list += generate_combinations(current_group)

                current_group['group_id'] = group_id
                current_group['cele'] = []
                current_group['hsap'] = []

            if taxonomy_id == 6239:
                current_group['cele'].append(gene_info)
            elif taxonomy_id == 9606:
                current_group['hsap'].append(gene_info)

        if (len(current_group['cele']) > 0) and (len(current_group['hsap']) > 0):
            ortholog_list += generate_combinations(current_group)

    with open(ortholog_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['entrez_id_cele', 'entrez_id_hsap'])
        for row in ortholog_list:
            writer.writerow(row)

def _get_homologene_entrez_wb_map(homologene_entrez_list_ce):
    """Returns the Entrez to WB ID map for Homologene.

    The first portion is obtained by using gene info table from
    <ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz>,
    downloaded on 2017-01-28. The ones that could not be found here are further found by
    scraping the ID history pages.

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between Entrez and WB IDs.
    """
    entrez_wb_df = pd.read_csv('data/entrez/Caenorhabditis_elegans.gene_info.gz', sep='\t',
                               header=0, usecols=[1, 5], names=['CE_ENTREZ', 'CE_WB_OLD'])

    # Pick out WB ID entries and separate with commas
    entrez_wb_df['CE_WB_OLD'] = entrez_wb_df['CE_WB_OLD'].apply(_get_wb_id_for_entrez)

    # Split comma-separated values in one row into multiple rows
    entrez_wb_df = tidy_split(entrez_wb_df, 'CE_WB_OLD', sep=',')

    # Find out the missing ones to scrape
    to_scrape = set(homologene_entrez_list_ce).difference(set(entrez_wb_df['CE_ENTREZ']))
    print("Number of genes to scrape: {}".format(len(to_scrape)))
    pd.Series(list(to_scrape)).sort_values() \
        .to_csv('data/homologene/entrez_wb_map_to_scrape.csv', index=False)

    # Load the rest from the scraped results
    scraped_df = pd.read_csv('data/homologene/entrez_wb_map_scraped.csv', sep=',',
                             names=['CE_ENTREZ', 'CE_WB_OLD'])

    entrez_wb_df = pd.concat([entrez_wb_df, scraped_df], axis=0)

    return entrez_wb_df

def _get_homologene_entrez_ensembl_map(homologene_entrez_list_hs):
    """Returns the Entrez to Ensembl map for Homologene.

    The first portion is obtained by using gene info table from
    <ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz>,
    downloaded on 2017-01-28. The ones that could not be found here are further found by
    scraping the ID history pages.

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between Entrez and Ensembl IDs.
    """
    entrez_ensg_df = pd.read_csv('data/entrez/Homo_sapiens.gene_info.gz', sep='\t',
                                 header=0, usecols=[1, 5], names=['HS_ENTREZ', 'HS_ENSG'])

    # Pick out WB ID entries and separate with commas
    entrez_ensg_df['HS_ENSG'] = entrez_ensg_df['HS_ENSG'].apply(_get_ensg_id_for_entrez)

    # Split comma-separated values in one row into multiple rows
    entrez_ensg_df = tidy_split(entrez_ensg_df, 'HS_ENSG', sep=',')

    # Find out the missing ones to scrape
    to_scrape = set(homologene_entrez_list_hs).difference(set(entrez_ensg_df['HS_ENTREZ']))
    print("Number of genes to scrape: {}".format(len(to_scrape)))
    pd.Series(list(to_scrape)).sort_values() \
        .to_csv('data/homologene/entrez_ensembl_map_to_scrape.csv', index=False)

    # Load the rest from the scraped results
    scraped_df = pd.read_csv('data/homologene/entrez_ensembl_map_scraped.csv', sep=',',
                             names=['HS_ENTREZ', 'HS_ENSG'])

    entrez_ensg_df = pd.concat([entrez_ensg_df, scraped_df], axis=0)

    return entrez_ensg_df

def _get_wb_id_for_entrez(entry):
    """Finds all the WB ID entries and returns them separated with commas."""
    id_set = re.findall(r'WormBase:(WBGene[0-9]{8})', entry)
    if len(id_set):
        return ','.join(id_set)
    else:
        return None

def _get_ensg_id_for_entrez(entry):
    """Finds all the ENSG entries and returns them separated with commas."""
    id_set = re.findall(r'ENSG[0-9]{11}', entry)
    if len(id_set):
        return ','.join(id_set)
    else:
        return None

####################
# WormBase
####################
def get_wormbase_table():
    """Table for WormBase ID to common name mapping.

    WormBase gene ID and InterPro domain tables (2016-09-15) from:
    <ftp://ftp.wormbase.org/pub/wormbase/releases/WS255/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz>
    <ftp://ftp.wormbase.org/pub/wormbase/releases/WS255/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS255.protein_domains.tsv>

    This table contains WormBase IDs, common names, as well as locus ID information.
    Ahringer RNAi clone locations (WS239) are read, updated, and left-joined
    with the WormBase table to provide non-ortholog C. elegans information.

    Returns:
        pandas.DataFrame: DataFrame containing mapping information between WB IDs, common names.
            locus IDs, and Ahringer RNAi clone locations.
    """
    wb_df = pd.read_csv("data/wormbase/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz", sep=',',
                        compression='gzip', header=None, usecols=[1, 2, 3],
                        names=['CE_WB_CURRENT', 'COMMON_NAME', 'LOCUS_ID'])

    # Read Ahringer locations, mapped to WS239
    ahringer_df = pd.read_csv("data/ahringer/locations_ws239.csv", sep=',',
                              header=None, usecols=[0,1],
                              names=["AHRINGER_LOC", 'CE_WB_OLD'])

    # Deal with WB ID changes
    ahringer_df = pd.concat([ahringer_df, get_ce_wb_updated(ahringer_df)], axis=1) \
            .sort_values(['CE_WB_CURRENT', 'AHRINGER_LOC']).reset_index(drop=True)

    # Read InterPro domains
    interpro_tuples = []
    with gzip.open('data/wormbase/c_elegans.PRJNA13758.WS255.protein_domains.tsv.gz', 'rt') as file:
        reader = csv.reader(file, delimiter="\t")
        for line in reader:
            wb_id = line[0]
            domains = line[3:]

            # InterPro domains are separated by '|' for each WormBase ID
            if len(domains):
                interpro_tuples.append((wb_id, "|".join(domains)))

    interpro_df = pd.DataFrame(interpro_tuples, columns=['CE_WB_CURRENT', 'INTERPRO_DOM'])

    # Left join using the WormBase gene ID table
    wb_df = pd.merge(wb_df, ahringer_df, how='left', on='CE_WB_CURRENT')
    wb_df = pd.merge(wb_df, interpro_df, how='left', on='CE_WB_CURRENT')

    wb_df = wb_df[['CE_WB_CURRENT', 'COMMON_NAME', 'LOCUS_ID', 'AHRINGER_LOC', 'INTERPRO_DOM']]

    return wb_df

if __name__ == "__main__":
    # Process the orthologs
    print('\nProcessing the orthologs...')
    ORTHOMCL = get_orthomcl()
    OMA = get_oma()
    COMPARA = get_compara()
    INPARANOID = get_inparanoid()
    ORTHOINSPECTOR = get_orthoinspector(preprocessed=True, uniprot_map_built=False)
    HOMOLOGENE = get_homologene(preprocessed=False)
    WORMBASE = get_wormbase_table()

    # Write to CSV
    print('\nWriting to CSV...')
    ORTHOMCL.to_csv('results/orthomcl.csv', index=False)
    OMA.to_csv('results/oma.csv', index=False)
    COMPARA.to_csv('results/compara.csv', index=False)
    INPARANOID.to_csv('results/inparanoid.csv', index=False)
    ORTHOINSPECTOR.to_csv('results/orthoinspector.csv', index=False)
    HOMOLOGENE.to_csv('results/homologene.csv', index=False)
    WORMBASE.to_csv('results/wormbase.csv', index=False)
    print('Done!')

    # Write to Excel
    print("\nWriting to Excel...")
    WRITER = pd.ExcelWriter('results/results.xlsx')
    ORTHOMCL.to_excel(WRITER, 'orthomcl', index=False)
    OMA.to_excel(WRITER, 'oma', index=False)
    COMPARA.to_excel(WRITER, 'compara', index=False)
    INPARANOID.to_excel(WRITER, 'inparanoid', index=False)
    ORTHOINSPECTOR.to_excel(WRITER, 'orthoinspector', index=False)
    HOMOLOGENE.to_excel(WRITER, 'homologene', index=False)
    WORMBASE.to_excel(WRITER, 'wormbase', index=False)
    WRITER.save()
    print('Done!')

    # Connect to database, create if needed
    print('\nConnecting to database...')
    ENGINE = create_engine('mysql://root:@localhost/ortholist')
    if not database_exists(ENGINE.url):
        create_database(ENGINE.url)

    # Write to database
    print('Writing to database \'ortholist\'...')
    ORTHOMCL.to_sql(name='orthomcl', con=ENGINE, if_exists='replace', index=False)
    OMA.to_sql(name='oma', con=ENGINE, if_exists='replace', index=False)
    COMPARA.to_sql(name='compara', con=ENGINE, if_exists='replace', index=False)
    INPARANOID.to_sql(name='inparanoid', con=ENGINE, if_exists='replace', index=False)
    ORTHOINSPECTOR.to_sql(name='orthoinspector', con=ENGINE, if_exists='replace', index=False)
    HOMOLOGENE.to_sql(name='homologene', con=ENGINE, if_exists='replace', index=False)
    WORMBASE.to_sql(name='wormbase', con=ENGINE, if_exists='replace', index=False)
    print('Done!')
