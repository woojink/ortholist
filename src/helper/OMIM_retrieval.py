#!/usr/bin/env python3

import csv
import gzip
from collections import defaultdict
import pandas as pd

def get_OMIM_annotations():

	"""
	Retrieves OMIM information from Ensembl 89

	Returns: A DataFrame containing SMART, GO, and HGNC annotations from Ensembl 89

	"""

	ensembl_data = defaultdict(lambda: defaultdict(set))

	with gzip.open('data/mart_export.txt.gz', 'rt') as file:

		reader = csv.reader(file, delimiter="\t")

		for ensg, omim_acc, omim_desc in reader:

			if omim_acc:

				ensembl_data[ensg]['OMIM_Accession'].add(omim_acc)

			if omim_desc:

				ensembl_data[ensg]['OMIM_Description'].add(omim_desc)


	ensembl_df = pd.DataFrame.from_dict(ensembl_data, orient='index')

	ensembl_df.reset_index(inplace=True)

	return ensembl_df

