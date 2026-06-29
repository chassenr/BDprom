#!/usr/bin/env Rscript

# to identify the candidate genepairs from gtf file

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
silent <- suppressMessages(lapply(package.list, function(X) {
  require(X, character.only = TRUE)
}))
rm(silent)
cat(" done\n")

#####################################################################
### Reading command line options ####

# define command line options
option_list <- list(
  # make_options for the workflow 1 (identification)
  make_option(
    c("-s", "--source"),
    type = "character",
    default = "src",
    help = "The directory containing Rscript for creating functions to extract candidate genes",
    metavar = "character"
  ),
  make_option(
    c("-i", "--indir"),
    type = "character",
    default = ".",
    help = "The input directory containing the parsed files",
    metavar = "character"
  ),
  make_option(
    c("-g", "--gtf"),
    type = "character",
    default = NULL,
    help = "The gtf input file",
    metavar = "character"
  ),
  make_option(
    c("-d", "--dist_max"),
    type = "integer",
    default = 600,
    help = "The maximum exon start position distance between two genes",
    metavar = "integer"
  ),
  # make_options for the workflow 2 (confirmation)
  make_option(
    c("-o", "--ge_dir"),
    type = "character",
    default = ".",
    help = "The directory containing the gene expression data",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cor_min"),
    type = "integer",
    default = 0.5,
    help = "The minimum correlation coeficcient between genes",
    metavar = "integer"
  ),
  make_option(
    c("-l", "--slope_diff"),
    type = "integer",
    default = 2,
    help = "The maximum fold deviation from a slope of 1 in the linear regression",
    metavar = "integer"
  ),
  make_option(
    c("-m", "--median_diff"),
    type = "integer",
    default = 1.2,
    help = "The maximum difference between median gene expression relative to the average of the medians",
    metavar = "integer"
  ),
  make_option(
    c("-t", "--detect_f "),
    type = "integer",
    default = 0.5,
    help = "The minimum fraction of studies the gene pair passes the above two conditions",
    metavar = "integer"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# confirm existence of makeoptions of workflow 1 and 2
if (is.null(opt$source) | is.null(opt$indir) | is.null(opt$gtf) | is.null(opt$dist_max) | is.null(opt$ge_dir) | is.null(opt$cor_min) | is.null(opt$slope_diff) | is.null(opt$median_diff) | is.null(opt$detect_f)) {
  print_help(opt_parser)
  stop("You must provide the input directory containing parsed files, the path to the Rscript, the gtf file, maximum distance between exons of two genes, and the name of the output file.\n", call. = FALSE)
}

# Set the working directory
# this will be taken out once script is made executable
# setwd("P:/CQ_Beratug_Bildung_GmbH/Internship_2026/IOW_2026/BDP_IOW/data")

# load the workspace
# load(identity_pairs.RData)
# check the input files
# file_list_hst2 <- list.files(opt$indir, pattern = "")

## Run workflow 1
# Set the source of the R script
source(paste0(opt$source, "/identify_candidate_gene_pairs.R"))

# Call the function
list_candidate_pairs <- identify_candidate_gene_pairs(
  gtf = opt$gtf,
  dist_max = opt$dist_max,
  out_file = paste0(opt$out_file, "/bdprom_candidate_pairs.txt")
)


## project part 3 #confirm the candidate genepairs
# we use the same source(R scripts) and output make options as the previous above code, since we have differnet file names

## Run workflow 2
# Set the source of the R script
source(paste0(opt$source, "/src/confirm_candidate_gene_pairs.R"))

# Call the function
ge_candidate_pairs <- confirm_candidate_gene_pairs(
  candidate_pairs = list_candidate_pairs,
  prior_out = opt$priorout,
  ge_dir = opt$gedir,
  cor_min = opt$cor,
  slope_diff = opt$slope,
  median_diff = opt$median,
  detect_f = opt$detect,
  out_file = paste0(opt$out_file, "/bdprom_ge_pairs.txt")
)

## Visualization
# Set the source of the R script
source(paste0("src/visualize_candidate_gene_pairs.R"))

# Call the function
visualize_candidate_gene_pairs(
  ge = ge_candidate_pairs,
  col_pal = rainbow,
  out_file = paste0(opt$out_file, "/bdprom_plot_pairs.pdf")
)

# save workspace
# save.image(identity_pairs.RData)
