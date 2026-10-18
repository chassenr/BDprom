#!/usr/bin/env Rscript


### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "tidyverse",
  "reshape"
)

# Function to check if packages are installed
is.installed <- function(pkg) {
  is.element(pkg, installed.packages()[, 1])
}

# If not all packages are available
if (any(!is.installed(package.list))) {
  cat("Not all required packages are available. They will now be installed.\n")
  # give the user the chance to abort manually
  Sys.sleep(20)

  # then install packages
  for (i in which(!is.installed(package.list))) {
    suppressMessages(install.packages(package.list[i], repos = "http://cran.us.r-project.org"))
  }
}

# Break the script if the package installation was unsuccessful
if (any(!is.installed(package.list))) {
  cat(
    paste0(
      "Unable to install all required packages.\nPlease install ",
      paste0(package.list[!is.installed(package.list)], collapse = ", "),
      " manually."
    )
  )
  break
}

# Load packages
cat("Loading libraries...")
silent <- suppressMessages(lapply(package.list, function(X) { require(X, character.only = TRUE) }))
rm(silent)
cat(" done\n")

# Some functions for message output
msg <- function(X) {
  cat(crayon::white(paste0("[", format(Sys.time(), "%T"), "]")), X)
}


### Reading command line options ####

# define command line options
option_list <- list(
  make_option(
    c("-s", "--sdir"),
    type = "character",
    default = "src",
    help = "The directory containing Rscript for creating functions to extract candidate genes",
    metavar = "character"
  ),
  make_option(
    c("-w", "--wdir"),
    type = "character",
    default = ".",
    help = "The input directory containing the gene expression count files",
    metavar = "character"
  ),
  make_option(
    c("-d", "--dist_max"),
    type = "integer",
    default = 600,
    help = "The maximum exon start position distance between two genes",
    metavar = "integer"
  ),
  make_option(
    c("-c", "--cor_min"),
    type = "numeric",
    default = 0.5,
    help = "The minimum pearson correlation coeficcient between genes",
    metavar = "numeric"
  ),
  make_option(
    c("-l", "--slope_diff"),
    type = "numeric",
    default = 2,
    help = "The maximum fold deviation from a slope of 1 in the linear regression",
    metavar = "numeric"
  ),
  make_option(
    c("-m", "--median_diff"),
    type = "numeric",
    default = 1.2,
    help = "The maximum difference between median gene expression relative to the average of the medians",
    metavar = "numeric"
  ),
  make_option(
    c("-f", "--detect_f"),
    type = "numeric",
    default = 0.5,
    help = "The minimum fraction of studies the gene pair passes the above two conditions",
    metavar = "numeric"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "Pairs_out",
    help = "Output directory name",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# confirm existence of makeoptions of workflow 1 and 2
# this is not necessary here, since script will work with defaults


### Load functions ####

source(paste0(opt$sdir, "/identify_candidate_gene_pairs.R"))
source(paste0(opt$sdir, "/confirm_candidate_gene_pairs.R"))
source(paste0(opt$sdir, "/visualize_candidate_gene_pairs.R"))


### Running the workflow ####

# Identification of candidate gene pairs
list_candidate_pairs <- identify_candidate_gene_pairs(
  gtf = paste0(opt$wdir, "/Genome/genome.gtf"),
  dist_max = opt$dist_max,
  out_file = paste0(opt$outdir, "/bdprom_candidate_pairs.txt")
)

# Confirmation of co-regulation based on gene expression
ge_candidate_pairs <- confirm_candidate_gene_pairs(
  candidate_pairs = list_candidate_pairs,
  prior_out = NULL,
  ge_dir = paste0(opt$wdir, "/Parsed_counts"),
  cor_min = opt$cor_min,
  slope_diff = opt$slope_diff,
  median_diff = opt$median_diff,
  detect_f = opt$detect_f,
  out_file = paste0(opt$outdir, "/bdprom_ge_pairs.txt")
)

# Visualization of co-expression
visualize_candidate_gene_pairs(
  ge = ge_candidate_pairs,
  col_pal = rainbow,
  out_file = paste0(opt$outdir, "/bdprom_plot_pairs.pdf")
)

# save workspace
save.image(paste0(opt$outdir, "/identify_wf.RData"))


