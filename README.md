# BDprom
Identification of bidirectional promoters in eukaryotic genomes based on gene expression profiles.

The approach consist of 3 workflows that can be run independently:

1. ```reads_wf```: Identification, download and processing of transcriptomic read data to obtain gene expression profiles
2. ```identify_wf```: Identification of candidate gene pairs likely controlled by a bidirectional promoter based on co-expression patterns
3. ```extract_wf```: Extraction, alignment, and motif identification of putative promoter sequences


# 1. reads_wf

This workflow automates the selection and processing of suitable RNA-seq data for an organism of your choosing from the Sequence Read Archive (SRA) via the ENA portal API to obtain gene expression profiles.

### Overview

<img src="https://github.com/chassenr/BDprom/blob/main/images/schematic_reads_wf.png" width="80%" height="80%">

**Obtaining metadata**: Taking the ```tax_id``` from the config file(see [usage](https://github.com/chassenr/BDprom/tree/main#usage)), this step automatically searches the SRA for RNA-seq experiments using the ENA portal API and downloads their metadata. Next, RNA-seq experiments are filtered to only retain those from comparable sequencing technology, with sufficient sequencing depth, replication of experimental conditions, and (optionally) those published in a peer-reviewed scientific article. The filter settings are user-configurable.

**Downloading fastq files**: Sequence read data is downloaded for each run from the ENA ftp site using [aria2c](https://aria2.github.io/manual/en/html/aria2c.html). After download MD5sum are compared and the workflow will abort if the checksums don't match. In those cases, the affected files will have to be manually removed and the workflow restarted. Successful downloads will not be repeated.

**Read quality control**: Sequence read data is quality cleaned by removing phiX reads, clipping sequencing adapters, and quality trimming with [BBduk](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/). The minimum read length for quality trimming is automatically determined based on the read length. For diagnostic purposes [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) profiles are generated for each run.

**Read mapping**: The reference genome for the organism of interest is downloaded and indexed for mapping with [Hisat2](https://daehwankimlab.github.io/hisat2/) in an exon- and splice-aware format. Sequence reads are mapped and bam files are filtered to retain only primary alignments and (in case of paired-end data) those mapping as proper pair.

**Gene expression quantification**: Reads (fragments for paired-end data) mapping to features based on the gene coordinated in the GTF file are counted with [FeatureCounts](https://subread.sourceforge.net/featureCounts.html) and aggregated in comma-separated values tables per study.

**Summarizing log information**: For each analysis step (phiX removal, adapter clipping, quality trimming, mapping, feature counting) the number and percentage of reads is summarized across all runs. It is highly recommended to review these tables to assess workflow performance.


### Usage

**Software installation**

This workflow uses [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). It has been developed and tested with snakemake version 8.28.0. Each rule requires a dedicated conda environment. At the moment, these conda environments are not automatically installed by snakemake, but assumed to the be available on your system. However, the yaml files for each environment can be found in the ```envs``` subdirectory of this repository. You further need ```wget```, [```dos2unix```](https://packages.debian.org/search?keywords=dos2unix) and [GNU ```parallel```](https://www.gnu.org/software/parallel/man.html) available on your system.

To get started, clone this repository. It is highly recommended to work with [SSH keys for cloning](https://docs.github.com/en/authentication/connecting-to-github-with-ssh).

```
git clone git@github.com:chassenr/BDprom.git
cd BDprom
```

**Workflow configuration**

The subdirectory ```config``` contains configuration files for each of the 3 workflows within this project. For the ```reads_wf``` you need to modify [```config_reads_wf.yaml```](https://github.com/chassenr/BDprom/blob/main/config/config_reads_wf.yaml). Example configuration settings for several species are available in ```example_config```().

* **```wdir```**: Path to your working directory. This is not the location of your repository clone, but the location where the output from the workflow will be written. You need at least **5 times the disk space** available than written in the log file for the ENA metadata retrieval and filtering (located at ```wdir/ENA_out/logs/filter_ena_metadata.log```)
* **```threads```** & **```parallel_threads```**: Number of default threads (should always be 1) and number of threads for parallelized analysis steps. Keep in mind that more is not always better. Ideally, the number of threads available on your system can be divided by ```parallel_threads```, so that multiple samples can be processed in parallel, each using ```parallel_threads```. A number between 16 and 48 should perform well. 
* **```gtf_link```**: Download link for the GTF genome annotation file. E.g. ```"http://ftp.ensemblgenomes.org/pub/protists/current/gtf/thalassiosira_pseudonana/Thalassiosira_pseudonana.ASM14940v2.62.gtf.gz"```
* **```genome_link```**: Download link for the genome fasta file. E.g. ```"http://ftp.ensemblgenomes.org/pub/protists/current/fasta/thalassiosira_pseudonana/dna/Thalassiosira_pseudonana.ASM14940v2.dna.toplevel.fa.gz"```. When using Ensembl resources, you should use the toplevel genome file without masking.
* **```tax_id```**: [NCBI taxon ID](https://www.ncbi.nlm.nih.gov/taxonomy) corresponding to the taxon for which the genome and gtf download links were specified. E.g. ```"35128"``` for Thalassiosira pseudonana.
* **```min_read_counts```**: Minimum number of reads per run. Runs with fewer reads will not be used in the analysis. E.g. 5000000.
* **```seq_technology```**: Allowed values for [sequencing platform](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#permitted-values-for-platform) (comma-separated list, no spaces). Only sequencing platforms generating comparable data should be listed here. E.g. "ILLUMINA,BGISEQ,DNBSEQ".
* **```min_n_runs```**: Minium number of sequencing runs per study. Studies with fewer runs after applying previous 2 filtering steps will not be used in the analysis. E.g. 5.
* **```only_keep_pubmed```**: Logical flag specifying if only ENA studies which can be linked to a scientific article on [PubMed](https://pubmed.ncbi.nlm.nih.gov/) or [PMC](https://pmc.ncbi.nlm.nih.gov/) should be analysed. Allowed values: "FALSE" or "TRUE".


**Workflow execution**

```
conda activate snakemake-8.28.0
snakemake -s rules/rules_reads_wf.smk --software-deployment-method conda --cores 384 -np
```

It is highly recommended to conduct a dry-run first (```-np```). Remove those flags to start the workflow.


### Output

The main output files of this worflow are:

* ```ENA_out/RNAseq_filtered.tsv```: A tab-separated values table with the metadata of the RNAseq runs, including PMC and PubMed accession number of related articles.
* ```Parsed_counts/gene_expression_*```: Comma-separated gene expression count tables, 1 file per ENA study. The first column contains the ```Ensembl_gene_id```, all subsequent collumns the read counts per feature and run.
* ```Abstracts/pmc_abstracts.txt``` & ```Abstracts/pubmed_abstracts.txt```: Titles and abstracts of the scientific articles linked to the RNAseq data in a tab-separated values table. This output is optional depending on the setting of ```only_keep_pubmed``` in the config file. 
* ```Logs_summary/summary_*```: Diagnostic information, summarizing the number and percentage of reads retained after each analysis step.


### Future development

* add conda env installation to snakemake rules
* add metadata filter based on abstract content


# 2. identify_wf

This workflow uses gene coordinates and co-expression patterns to create a list of gene pairs potentially regulated by a bidirectional promoter. It further provides a visual summary of the co-expression for manual refinement of the list of candidate gene pairs.

So far, the ```identify_wf``` workflow is only available as interactive R script using the commands outlined below in this README.

### Preparations

**Obtain genome coordinates and differential gene expression data:**

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

If gene expression data was obtained via the ```reads_wf``` workflow, the files ```gene_expression_*``` in the ```Parsed_counts``` output directory will be in the correct format and can be used directly. Any other source of gene expression data will also be accepted as long as file name and content meet the above criteria.


**Prepare R working environment:**

Apart from the genome coordinates and gene expression data, you will need the R functions provided in this repository. For this, clone or download the repository to your local computer. The below instructions assume that the gtf file and the directory with the gene expression tables are located in the main repository directory, however, you can easily provide file paths to other locations. 

```
setwd("path/to/your/working/directory")
require(tidyverse)
require(reshape)
```

### Identification of candidate gene pairs based on genome organisation

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


### Validation of co-regulation based on gene expression profiles 

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


### Visualization of co-expression

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


### Future development

* There will be an automated version of the ```identify_wf``` workflow available soon, implemented in snakemake.


# 3. extract_wf

:warning: WIP :warning:


# Citation

* WIP

Please also cite the dependencies:

* ENA portal API
* R Core Team (2024) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.
* Davis T (2024) optparse: Command Line Option Parser. R package version 1.7.5. https://CRAN.R-project.org/package=optparse.
* Allaire J (2023) config: Manage Environment Specific Configuration Values. R package version 0.3.2. https://CRAN.R-project.org/package=config.
* Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019) Welcome to the tidyverse. Journal of Open Source Software 4(43):1686. https://doi.org/10.21105/joss.01686.
* Bengtsson H (2023) R.utils: Various Programming Utilities. R package version 2.12.3. https://CRAN.R-project.org/package=R.utils.
* NCBI edirect
* aria2
* Bushnell B (2014) BBMap. Available online at https://sourceforge.net/projects/bbmap/.
* Andrews S, Krueger F, Segonds-Pichon A, Biggins L, Krueger C, Wingett S (2010) FastQC: a quality control tool for high throughput sequence data. Available online at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
* Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019) Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 37(8):907-915. https://www.nature.com/articles/s41587-019-0201-4.
* Liao Y, Smyth GK, Shi W (2014) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30(7):923-30. https://doi.org/10.1093/bioinformatics/btt656.
