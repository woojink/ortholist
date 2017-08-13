import csv
import gzip

from collections import defaultdict

import pandas as pd

from databases.Database import Database
from helper.misc import generate_combinations

class OrthoMCL(Database):
    """OrthoMCL ortholog table between C. elegans and H. sapiens.

    The H. sapiens side is provided as ENSP IDs instead of ENSG IDs, so the
    IDs are mapped using Ensembl release 56, which OrthoMCL version 5 uses.

    Version 5 (2015-07-23)
    """
    def __init__(self):
        super().__init__(name="OrthoMCL", filename="orthomcl")

    def _read_raw(self):
        """Returns an ortholog table for OrthoMCL

        Because OrthoMCL orthologs are provided as groupings, combinations are
        generated using generate_combinations(), which uses itertools.product().

        Returns:
            A DataFrame containing the raw orthologs from OrthoMCL
        """
        ortholog_list = []
        with gzip.open('data/orthomcl/groupings.csv.gz', 'rt') as file:
            reader = csv.reader(file, delimiter=',')

            groups = defaultdict(lambda: defaultdict(set))
            for row in reader:
                group_id = row[1]
                source_id = row[0]

                if source_id.startswith('WBGene'):
                    groups[group_id]['cele'].add(source_id)
                else:
                    groups[group_id]['hsap'].add(source_id)

        for key in groups:
            ortholog_list += generate_combinations(groups[key])

        return pd.DataFrame(ortholog_list, columns=['CE_WB_OLD', 'HS_ENSP']) \
                .drop_duplicates()

    def _perform_worm_mapping(self):
        return pd.merge(self.df, self._get_ensembl_56_ensp_ensg_map(),
                        left_on='HS_ENSP', right_index=True)

    @staticmethod
    def _get_ensembl_56_ensp_ensg_map():
        """Returns ENSP to ENSG mapping from Ensembl release 56

        Because there convinient direct mapping from ENSP to ENSG, various
        tables are used to map between the different stages.

        Data from
            ftp://ftp.ensembl.org/pub/release-56/mysql/homo_sapiens_core_56_37a/

        The schema is described in
            data/ensembl/56/homo_sapiens_core_56_37a.sql.gz

        Returns:
            A DataFrame containing mapping between ENSP and ENSG
        """
        ensp_ensg_df = pd.read_csv( \
                        "data/ensembl/56/translation_stable_id.txt.gz",
                        sep='\t', compression='gzip',
                        header=None, usecols=[0, 1],
                        names=["translation_id", "ENSP"])

        translation_to_transcript = pd.read_csv( \
                        "data/ensembl/56/translation.txt.gz", sep='\t',
                        compression='gzip', header=None, usecols=[0, 1],
                        names=["translation_id", "transcript_id"])
        transcript_to_gene = pd.read_csv( \
                        "data/ensembl/56/transcript.txt.gz", sep='\t',
                        compression='gzip', header=None, usecols=[0, 1],
                        names=["transcript_id", "gene_id"])
        gene_id_to_gene = pd.read_csv( \
                        "data/ensembl/56/gene_stable_id.txt.gz", sep='\t',
                        compression='gzip', header=None, usecols=[0, 1],
                        names=["gene_id", 'HS_ENSG'])

        ensp_ensg_df = ensp_ensg_df \
            .set_index("translation_id") \
            .join(translation_to_transcript \
            .set_index("translation_id")) \
            .set_index("transcript_id") \
            .join(transcript_to_gene.set_index("transcript_id")) \
            .set_index("gene_id") \
            .join(gene_id_to_gene.set_index("gene_id")) \
            .set_index("ENSP")

        return ensp_ensg_df
