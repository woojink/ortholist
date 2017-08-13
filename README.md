# Ortholist 2
## Overview
This project updates Ortholist ([version 1 here](http://www.greenwaldlab.org/ortholist/)) with results from both additional ortholog prediction programs and updated releases of previously included programs.

## Databases used
* Ensembl Compara (Releases 87, 2016-12; 88, 2017-03; and 89, 2017-05)
* Homologene (Build 68, 2014-05-06)
* InParanoid (8.0, 2013-12)
* OMA (May 2016)
* OrthoInspector ("Quest for Orthologs (QfO) 04_2013", 2013-04)
* OrthoMCL (Version 5, 2015-07-23)

## Instructions
### Pre-requisites
This project makes use of Python 3. If you're on macOS, you can get these by first installing [Homebrew](https://brew.sh/) by putting the following on your Terminal:
```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Then running
```bash
brew install python3
```

## Methodology
### WormBase
The different databases used for Ortholist were generated at different dates. As a result, worm genes are sometimes not up to date with the most current WormBase database. We keep track of every changed WormBase ID and update it to the most current version or to none if deprecated. The WormBase version used for this version of Ortholist is WS255, available on <`ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.WS255.geneIDs.txt.gz)`>.

This WormBase ID updating is done through [`get_ce_wb_updated()`](https://github.com/woojink/ortholist/blob/master/src/helper/wb_map.py#L49) under `helper.wb_map`

### Ensembl Compara
Releases 87 (2016-12), 88 (2017-03), and 89 (2017-05) of Ensembl are used for Ortholist.

1. Orthologs are directly obtained from [December 2016](http://dec2016.archive.ensembl.org/biomart/martview), [March 2017](http://mar2017.archive.ensembl.org/biomart/martview), and [May 2017 BioMart](http://may2017.archive.ensembl.org/biomart/martview). Worm genes were obtained as Wormbase IDs and human genes as ENSG IDs.
2. WormBase ID changes are dealt with using `get_ce_wb_updated()` (see [above](#wormbase))

#### Statistics (combined 87-89)
* Before
    - Ortholog pairs: 16,340
    - Unique worm genes: 6,825
    - Unique human genes: 9,188
* After
    - Orthologs pairs: 16,339
    - Unique worm genes: 6,823 (1 merged, 1 pseudogene entries)
    - Unique human genes: 9,188

### Homologene
Build 68 (2014-05-06) is used for Homologene

1. Raw SQL table (`homologene.data`) for the orthologs is downloaded from <`ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data`>. Both worm and human genes are provided as Entrez IDs. As described on Homologene's README located on <`ftp://ftp.ncbi.nih.gov/pub/HomoloGene/README`>, the first three columns of the provided database is as follows:

    1. HID (HomoloGene group id)
    2. Taxonomy ID
    3. Gene ID

    Taxonomy IDs for _C. elegans_ and _H. sapiens_ are 6239 and 9606, respectively. We first filter out any ortholog that isn't worm or human with [`preprocess.sh`](https://github.com/woojink/ortholist/blob/master/data/homologene/preprocess.sh).

2. The orthologs are provided as groupings similar to InParanoid, so combinations are generated using [`generate_combinations()`](https://github.com/woojink/ortholist/blob/master/src/helper/misc.py#L41) under `helper.misc`. We first generate a more tractable ortholog file `homologene.tsv` under the Homologene data folder using the above method.
3. Entrez IDs need to be converted to WormBase and Ensembl IDs. We first obtain the Entrez gene info table from <`ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz`> (2017-01-28). There are some entries that correspond to multiple WormBase IDs, which are separated into WormBase entities (see [`tidy_split()`](https://github.com/woojink/ortholist/blob/master/src/helper/misc.py) under `helper.misc` for the row splitting method).
4. There are entries now missing in Homologene since 2014-05-06, so we scrape the history pages as we had done for UniProt entries. The full history page URLs are provided with the following format:
    ```
    https://www.ncbi.nlm.nih.gov/gene/{Entrez ID}?report=xml&format=text
    ```

    WormBase IDs are extracted from these pages. The remaining IDs that could not be determined using this method are obtained through manually searching WormBase itself. The results from all three methods allow Entrez to WormBase ID mapping.

5. Likewise with human genes, the Entrez gene info table is downloaded from <`ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz`> (2017-01-28) first. Scraping was attempted, but none of the entries returned an Ensembl ID, so the remaining 6 IDs were manually matched to Ensembl entries using names and other properties. These two methods provided Entrez to Ensembl ID mapping.
6. Lastly, the WormBase ID changes are dealt with using `get_ce_wb_updated()` (see [above](#wormbase))

#### Statistics
* Before
    - Orthologs pairs: 3,811
    - Unique worm genes: 3,799
    - Unique human genes: 3,184
* Entrez - Wormbase mapping (3,799 of 3,799)
    - Found through Entrez gene info table: 3,765
    - Found through scraping history pages: 24
    - Manually matched: 10
    - Unique Wormbase IDs: 3,785
* Entrez - Ensembl mapping (3,184 of 3,184)
    - Found through the ID mapping tool: 3,178
    - Found through scraping history pages: 0
    - Manually matched: 6
    - Unique Ensembl IDs: 3,205
* After mapping and Wormbase ID update
    - Orthologs pairs: 3,817
    - Unique worm genes: 3,779 (2 merged, 6 pseudogene entries)
    - Unique human genes: 3,205

### InParanoid
Release 8.0 (2013-12) of InParanoid is used for Ortholist.

1. Raw SQL table (`sqltable.C.elegans-H.sapiens`) for the orthologs is downloaded from [this link](http://inparanoid.sbc.su.se/download/8.0_current/Orthologs_other_formats/C.elegans/InParanoid.C.elegans-H.sapiens.tgz). Both worm and human genes are provided as UniProt IDs. The orthologs are provided as groupings, so combinations are generated using [`generate_combinations()`](https://github.com/woojink/ortholist/blob/master/src/helper/misc.py#L41) under `helper.misc`. We first generate a more tractable ortholog file `orthologs.tsv` under the InParanoid data folder using the above method.
2. For the worm genes:
    * The majority of the UniProt IDs are mapped using the ID mapping tool from UniProt available [here](http://www.uniprot.org/uploadlists/). As release 8.0 was released in 2013-12, there are many IDs no longer mappable.
    * For the rest, we look at the UniProt ID history pages to determine the last known entry (e.g. [`http://www.uniprot.org/uniprot/A4UVJ9?version=*`](http://www.uniprot.org/uniprot/A4UVJ9?version=*) shows 47 is the latest for `A4UVJ9`)
    * Then using that latest entry, we access the actual entry page to scrape the last known WormBase ID for the particular UniProt ID (e.g. [`http://www.uniprot.org/uniprot/A4UVJ9.txt?version=47`](http://www.uniprot.org/uniprot/A4UVJ9.txt?version=47) shows `WBGene00019439` is the last known correspondent with `A4UVJ9`)
3. For the human genes:
    * As with the worm genes, the majority of the UniProt IDs are mapped using the [ID mapping tool](http://www.uniprot.org/uploadlists/)
    * The ones that could not be found at this point are then searched through the [BioMart tool for Ensembl 74](http://dec2013.archive.ensembl.org/biomart/martview/) (2013-12), which is the closest release to InParanoid 8.0 (also 2013-12). Both SwissPort and TrEMBL entries are looked at.
    * At this point, the remaining 135 IDs are scraped through the history pages as with the worm genes previously, yield 9 more entries.
4. Lastly, the WormBase ID changes are dealt with using `get_ce_wb_updated()` (see [above](#wormbase))

#### Statistics
* Before
    - Orthologs pairs: 12,826
    - Unique worm genes: 5,545
    - Unique human genes: 8,395
* Uniprot - Wormbase mapping (5,545 of 5,545)
    - Found through the ID mapping tool: 5,482
    - Found through scraping history pages: 63
* Uniprot - Ensembl mapping (8,272 of 8,272)
    - Found through the ID mapping tool: 8,220
    - Found through the BioMart tool (SwissProt): 27
    - Found through the BioMart tool (TrEMBL): 13
    - Found through scraping history pages: 12
    - Remaining: 123
* After
    - Ortholog pairs: 15,136
    - Unique worm genes: 5,587
    - Unique human genes: 8,957

### OMA
May 2016 release of OMA is used.

1. Orthologs are directly obtained from [this endpoint](http://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=CAEEL&p2=HUMAN&p3=EnsemblGene). Each line of this file indicates an ortholog: worm genes are provided as OMA IDs and human genes as ENSG IDs.
2. OMA IDs for worm genes are first mapped to Wormpep IDs using the provided [map file](http://omabrowser.org/All/oma-wormbase.txt.gz) from OMA. The Wormpep IDs are from WS235 and they are translated to Wormbase IDs using the `wormpep.table235` from Wormbase <`ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/protein/c_elegans.WS235.wormpep_package.tar.gz`>.
3. WormBase ID changes are dealt with using `get_ce_wb_updated()` (see [above](#wormbase))

#### Statistics
* Before
    - Orthologs pairs: 5,370
    - Unique worm genes: 3,885
    - Unique human genes: 4,559
* Genes missing due to Wormpep to Wormbase mapping: 0
* After
    - Orthologs pairs: 5,368
    - Unique worm genes: 3,884 (9 merging and 2 pseudogene events)
    - Unique human genes: 4,559

### OrthoInspector
For OrthoInspector, "Quest for Orthologs (QfO) 04_2013" (2013-04) is used.

1. The raw worm orthologs are available on [this `gzip`-ed CSV file](http://lbgi.fr/orthoinspector/dbquery_qfo/data/CSV/62.csv.gz) (warning: large 1.22 GB file). Each line in this file represents an ortholog. Only the UniProt IDs are extracted from each line, then filtered for only lines containing human orthologs (see [`preprocess.sh`](https://github.com/woojink/ortholist/blob/master/data/orthoinspector/preprocess.sh) for more details)
2. The rest of the process for OrthoInspector is identical, where
    * Worm genes were first mapped through UniProt's ID mapping tool, then scraped through history pages
    * Human genes were also first mapped through UniProt's ID mapping tool, then the [BioMart tool for Ensembl 74](http://dec2013.archive.ensembl.org/biomart/martview/) (2013-12), which is the closest release to 2013-04, then lastly through scraping. 11 of 13 remaining ID correspondents are found after scraping.
3. Lastly, the WormBase ID changes are dealt with using `get_ce_wb_updated()` (see [above](#wormbase))

#### Statistics
* Before
    - Orthologs pairs: 10,542
    - Unique worm genes: 5,547
    - Unique human genes: 7,317
* Uniprot - Wormbase mapping (5,547 of 5,547)
    - Found through the ID mapping tool: 5,480
    - Found through scraping history pages: 67
    - Unique Wormbase IDs: 5,370
* Uniprot - Ensembl mapping (7,316 of 7,317)
    - Found through the ID mapping tool: 7,279
    - Found through the BioMart tool (SwissProt): 15
    - Found through the BioMart tool (TrEMBL): 10
    - Found through scraping history pages: 11
    - Manually matched: 1
    - Unmatched: 1 (pseudogene)
    - Unique Ensembl IDs: 7,782
* After mapping and Wormbase ID update
    - Orthologs pairs: 11,455
    - Unique worm genes: 5,364 (11 merged, 10 pseudogene entries)
    - Unique human genes: 7,782

### OrthoMCL
Version 5 (2015-07-23) of OrthoMCL is accessible [here](http://orthomcl.org/). Worm genes are provided as WormBase IDs and human genes are provided as ENSP IDs. Version 5 of OrthoMCL uses Ensembl release 56 <`ftp://ftp.ensembl.org/pub/release-56/mysql/homo_sapiens_core_56_37a/`>.

1. First, the orthologs are obtained from OrthoMCL with the following strategy:
    * Step 1: Phyletic filter with `cele+hsap=2T`
    * Step 2: To sequences
    * Step 3a: Nested strategy of union between `cele` and `hsap` taxonomies
    * Step 3b: Intersection between output of Step 2 and 3a
    * This results in 13,515 entries, with source accession IDs and group IDs
2. The orthologs are provided as groupings, so combinations are generated using [`generate_combinations()`](https://github.com/woojink/ortholist/blob/master/src/helper/misc.py#L41) under `helper.misc`.
3. Human genes are provided as ENSP IDs, which needs to be converted to ENSG IDs for the purpose of this project. Using four tables from Ensembl (`translation_stable_id.txt.gz`, `translation.txt.gz`, `transcript.txt.gz`, `gene_stable_id.txt.gz`), we are able to obtain ENSG ID that correspond to each ENSP ID.
4. WormBase ID changes are dealt with using `get_ce_wb_updated()` (see [above](#wormbase))

#### Statistics
* Before
    - Orthologs pairs: 12,718
    - Unique worm genes: 5,727
    - Unique human proteins: 7,788
* Genes missing due to `ENSP` to `ENSG` mapping: 0
* After
    - Orthologs pairs: 12,680
    - Unique worm genes: 5,707 (63 merged, 33 pseudogene, 4 killed due to lack of evidence, 1 transposon entries)
    - Unique human genes: 7,784


