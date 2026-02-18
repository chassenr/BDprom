# BDprom
Identification of bidirectional promoters in eukaryotic genomes based on gene expression profiles




## Step 0: Preparations

**Download genome coordinates and differential gene expression data:**

:warning: **WIP**: add download instructions

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
* exon start positions of the 2 genes in their respecitve reading direction are less than a pre-defined threshold (default: 600bp) distant from each other

As these criteria solely rely on genome organisation, only a gtf file is required for this analysis step. If no file name / path is provided to the command, R will look in the working directory for a matching file type.

:warning: **WIP**: show example of gtf file content.

The selection is then performed with the R function ```[identify_candidate_gene_pairs()](https://github.com/chassenr/BDprom/blob/main/src/identify_candidate_gene_pairs.R)```:

```
source("src/identify_candidate_gene_pairs.R")
list_candidate_pairs <- identify_candidate_gene_pairs(
  gtf = "Sadia_Phaeodactylum_tricornutum.ASM15095v2.59._agat.gtf",
  dist_max = 600,
  out_file = "bdprom_candidate_pairs.txt"
)
```

This command write the coordinates of the candidate gene pairs to file and creates a list in your workspace that contains both the full table of genes corrdinates for the organism as well as a subset thereof for the candidate gene pairs.

:warning: **WIP**: describe columns of output file


## Step 2: Validation of co-regulation based on gene expression profiles 

Second, gene expression profiles across multiple samples are used to confirm co-regulation of each gene pair. Ideally, data from multiple studies is used covering different experimental / environmental conditions to increase the chances of observing differential gene expression for all candidate pairs. Gene expression data can be obtained from [DiatomicBase](https://www.diatomicsbase.bio.ens.psl.eu/) as described above. For each study, read count data should be saved as csv. By default R will look for such files in a subdirectory ```diatomicbase_data``` in the current working directory. For the evaluation of co-expression patterns, TPM (transcripts per million) will then be calculated.

:warning: **WIP**: show example of csv read count file content.

The candidate gene pairs from step 1 are then filtered based on the following criteria:

* minimum correlation coefficient (default: 0.5)
* maximum fold deviation from a slope of 1 in the linear regression (default: 2, i.e. slope ranging between 0.5 and 2)
* minimum fraction of studies the gene pair passes the above two conditions (default: 0.5)

Gene expression patterns are then analysed with the R function ```[confirm_candidate_gene_pairs()](https://github.com/chassenr/BDprom/blob/main/src/confirm_candidate_gene_pairs.R)```:

```
source("src/confirm_candidate_gene_pairs.R")
ge_candidate_pairs <- confirm_candidate_gene_pairs(
  list_candidate_pairs,
  ge_dir = "diatomicbase_data",
  cor_min = 0.5
  slope_diff = 2,
  detect_f = 0.5
)
```

## Step 3: Visualization of co-expression

Additional to the summary statistics about the co-expression patterns generated in step 2, pairwise scatter plots are created for each gene pair passing the filter in step 2. For studies where the gene pair is expressed, but which did not make the thresholds in step 2, gene expression (TPM) is shown as grey dots. Otherwise points are colored by study and the correlation coefficient for the gene pair is included in the legend. Additionally, the line, slope and R2 of the linear regression using gene expression values from all studies (incl. the grey dots) and only using the values from the studies where the gene pair passed the filter are shown in grey and black, respectively.

The R plotting function tailored for this analysis is ```[visualize_candidate_gene_pairs()](https://github.com/chassenr/BDprom/blob/main/src/visualize_candidate_gene_pairs.R)```:

```
source("src/visualize_candidate_gene_pairs.R")
visualize_candidate_gene_pairs(
  ge_candidate_pairs,
  col_pal = rainbow,
  out_file = "bdprom_plot_pairs.pdf"
)
```
