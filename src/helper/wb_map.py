"""
This module provides helper function to provide current Wormbase ID mapping
for deprecated gene IDs. If the gene is not present in WS255, the current
status is read from Dan's mapping table.

WS255 list was downloaded from:
ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz
"""
import csv
import gzip

import pandas as pd

# Read genes from Wormbase WS255
WB_WS255 = set()
with gzip.open('data/wormbase/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz', \
               'rt') as f:
    READER = csv.reader(f, delimiter=",")
    for row in READER:
        WB_WS255.add(row[1])

# Read Dan's mapping for changed IDs
WB_OLD_TO_CURRENT_MAP = {}
with open('data/wormbase/WB_changes.csv') as f:
    READER = csv.reader(f, delimiter=",")
    next(READER)
    for old, current, comment in READER:
        if current == "":
            current = None
        WB_OLD_TO_CURRENT_MAP[old] = (current, comment)

def get_ce_wb_current(wb_id_old):
    """Provides the current WB ID if present"""
    if wb_id_old in WB_OLD_TO_CURRENT_MAP:
        return WB_OLD_TO_CURRENT_MAP[wb_id_old][0]
    elif wb_id_old in WB_WS255:
        return wb_id_old
    return None

def get_ce_wb_comment(wb_id_old):
    """Provides the comment for the changed IDs if present"""
    if wb_id_old in WB_OLD_TO_CURRENT_MAP:
        return WB_OLD_TO_CURRENT_MAP[wb_id_old][1]
    elif wb_id_old in WB_WS255:
        return None
    return "Not mapped, not in WS255"

def get_ce_wb_updated(df):
    """Returns a curated table with current IDs and comments

    Given a column of WB IDs, return a curated table with current IDs
    and the comment for the change (if applicable)
    """
    ids, comments = df['CE_WB_OLD'].apply(get_ce_wb_current), \
                        df['CE_WB_OLD'].apply(get_ce_wb_comment)
    new = pd.concat([ids, comments], axis=1)
    new.columns = ['CE_WB_CURRENT', 'CE_WB_COMMENT']

    return new
