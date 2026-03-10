# BDprom
Identification of bidirectional promoters in eukaryotic genomes based on gene expression profiles


## Step 0: Preparations

**Download genome coordinates and differential gene expression data:**

Gene annotations and coordinates can be obtained from [EnsemblGenomes](https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/current/gtf/) for the organisms available on this platform. E.g. Navigate to phaeodactylum_tricornutum, download Phaeodactylum_tricornutum.ASM15095v2.62.gtf.gz (e.g. into the working directory used for running the below commands), and unzip it.

Gene expression data is available on [DiatomicBase](https://www.diatomicsbase.bio.ens.psl.eu/). The below workflow expects that the file names for the gene expression count data end with the study accession number (PRJ*) immediately prior to the file extension and are saved as comma-separated values (.csv), e.g.: ```Converted_Counts_Data_PRJNA278661.csv```.

Each file with gene expression data must contain a column named "Ensembl_gene_id" with gene names. These gene names must correspond to the gene IDs in the gtf file. The columns in the csv file containing the gene expression count data are recognized automatically by selecting all columns with integer data. Please ensure that no integer columns other than from gene expression quantification are included.

```
"User_ID sorted by SD","Ensembl_gene_id","symbol","early_phase_rep1","early_phase_rep2","early_phase_rep3","stationary_phase_rep1","stationary_phase_rep2","stationary_phase_rep3"
" PHATR3_J41316","PHATR3_J41316","  ",3219,4170,3657,6,8,10
" PHATR3_J45944","PHATR3_J45944","  ",2990,4048,3469,15,17,16
" PHATR3_EG02608","PHATR3_EG02608","  ",1120,1460,1305,181963,197537,142801
" PHATR3_J46547","PHATR3_J46547","  ",7734,10183,9089,54,56,64
```


**Prepare R working environment:**

Apart from the genome coordinates and gene expression data, you will need the R functions provided in this repository. For this clone or download the repository to your local computer. The below instructions assume that the gtf file and the directory with the gene expression tables are located in the main repository directory, however, you can easily provide file paths to other locations. 

```
setwd("path/to/your/working/directory")
require(tidyverse)
require(reshape)
```

## Step 1: Identification of candidate gene pairs based on genome organisation

In a first step gene pairs, which may be regulated by one bidirectional promoter, are identified based on their position and orientation in the genome. In order to be selected the need to meet the following criteria:

* the 2 genes are read in opposite orientation
* they are neighbors of each other without any intervening genes
* they are located on the same chromosome
* exon start positions of the 2 genes in their respecitve reading direction are less than a pre-defined threshold (default: ```dist_max = 600```bp) distant from each other

As these criteria solely rely on genome organisation, only a gtf file is required for this analysis step (```gtf```). If no file name / path is provided to the command, R will look in the working directory for a matching file type.

The selection is then performed with the R function [```identify_candidate_gene_pairs()```](https://github.com/chassenr/BDprom/blob/main/src/identify_candidate_gene_pairs.R):

```
source("src/identify_candidate_gene_pairs.R")
list_candidate_pairs <- identify_candidate_gene_pairs(
  gtf = "Phaeodactylum_tricornutum.ASM15095v2.62.gtf",
  dist_max = 600,
  out_file = "bdprom_candidate_pairs.txt"
)
```

This command writes the coordinates of the candidate gene pairs to file (```out_file```) and creates a list in your workspace that contains both the full table of genes coordinates for the organism as well as a subset thereof for the candidate gene pairs.

```
(rowname)	length	min	max	orientation	start	chromosome	neighbor
PHATR3_J42424	1656	18350	20109	-	20109	1	PHATR3_J7430
PHATR3_J7430	1764	20334	22097	+	20334	1	PHATR3_J42424
PHATR3_EG02408	2700	24703	27789	-	27789	1	PHATR3_J31409
PHATR3_J31409	636	28204	28839	+	28204	1	PHATR3_EG02408
```

* length: length of the gene (cumulative exon length)
* min / max: first and last coordinates of gene (across all exons) in chromosome reading direction
* orientation: gene orientation relative to chromosome reading direction
* start: start position of first exon in gene in gene reading direction (i.e. start position of minus oriented genes will be higher than end position)
* chromosome: chromosome ID on which gene is located
* neighbor: second gene in candidate pair


## Step 2: Validation of co-regulation based on gene expression profiles 

Second, gene expression profiles across multiple samples are used to confirm co-regulation of each gene pair. Ideally, data from multiple studies is used covering different experimental / environmental conditions to increase the chances of observing differential gene expression for all candidate pairs. Gene expression data can be obtained from [DiatomicBase](https://www.diatomicsbase.bio.ens.psl.eu/) as described above. By default R will look for csv files with gene expression data in a subdirectory ```ge_dir = "diatomicbase_data"``` in the current working directory. For the evaluation of co-expression patterns, TPM (transcripts per million) will then be calculated.

The candidate gene pairs from step 1 (```candidate_pairs = list_candidate_pairs```) can then be filtered based on the following criteria:

* minimum correlation coefficient (default: ```cor_min = 0.5```)
* maximum fold deviation from a slope of 1 in the linear regression (default: ```slope_diff = 2```, i.e. slope ranging between 0.5 and 2)
* maximum difference between median gene expression relative to the average of the medians (default: ```median_diff = 1.5```)
* minimum fraction of studies the gene pair passes the above two conditions (default: ```detect_f = 0.5```)

It is possible to modify or disable any of these filters by changing the values in the below code. To disable a filter it should be set to ```NULL```. 

Gene expression patterns are then analysed with the R function [```confirm_candidate_gene_pairs()```](https://github.com/chassenr/BDprom/blob/main/src/confirm_candidate_gene_pairs.R):

```
source("src/confirm_candidate_gene_pairs.R")
ge_candidate_pairs <- confirm_candidate_gene_pairs(
  candidate_pairs = list_candidate_pairs,
  prior_out = NULL,
  ge_dir = "diatomicbase_data",
  cor_min = 0.5,
  slope_diff = 2,
  median_diff = 1.5,
  detect_f = 0.5,
  out_file = "bdprom_ge_pairs.txt"
)
```

The command writes the co-expression information for the gene pairs passing all applied filters to file (```out_file```) and creates a list in your workspace that additionally contains the gene expression data in TPM and the co-expression information for each candidate gene pair prior to filtering. Since the computation of co-expression patterns can be time-consuming, there is also the option to only repeat the filter without having to re-process the gene expression count tables by providing the output R object from the first run to ```prior_out```.

```
gene_pair	gene_1	gene_2	cor	slope	r2	min_1	min_2	median_1	median_2	max_1	max_2	study
PHATR3_EG00410___PHATR3_J48873	PHATR3_EG00410	PHATR3_J48873	0.962	0.975	0.971	0.083	0.063	2.984	2.541	6.670	7.290	PRJEB26173
PHATR3_EG00410___PHATR3_J48873	PHATR3_EG00410	PHATR3_J48873	0.810	0.603	0.978	17.781	13.973	25.654	14.786	31.0361	17.456	PRJEB34512
PHATR3_EG00410___PHATR3_J48873	PHATR3_EG00410	PHATR3_J48873	0.822	1.422	0.969	3.729	5.546	7.713	12.441	11.340	15.263	PRJNA278661
```

* gene_pair: identifer for the gene pair consisting of alphabetically sorted names of the two genes per pair separated by 3 underscores
* gene_1 / gene_2: names of the two genes per pair listed in separate columns
* cor: Pearson correlation coefficient between the TPM gene expression profiles of the two genes per pair
* slope: slope of the linear regression of gene_2 against gene_1 (intercept forced to zero)
* min_1 / min_2: minimum TPM in the data set (SRA study)
* median_1 / median_2: median TPM in the data set (SRA study)
* max_1 / max_2: maximum TPM in the data set (SRA study)
* study: SRA study accession for which co-expression was calculated


## Step 3: Visualization of co-expression

Additional to the summary statistics about the co-expression patterns generated in step 2, pairwise scatter plots are created for each gene pair passing the filter in step 2. For studies where the gene pair is expressed, but which did not make the thresholds in step 2, gene expression (TPM) is shown as grey dots. Otherwise points are colored by study and the correlation coefficient for the gene pair is included in the legend. Additionally, the line, slope and R2 of the linear regression using gene expression values from all studies (incl. the grey dots) and only using the values from the studies where the gene pair passed the filters are shown in grey and black, respectively.

The R plotting function tailored for this analysis is [```visualize_candidate_gene_pairs()```](https://github.com/chassenr/BDprom/blob/main/src/visualize_candidate_gene_pairs.R) using as input the R object created in step 2 (```ge = ge_candidate_pairs```):

```
source("src/visualize_candidate_gene_pairs.R")
visualize_candidate_gene_pairs(
  ge = ge_candidate_pairs,
  col_pal = rainbow,
  out_file = "bdprom_plot_pairs.pdf"
)
```
