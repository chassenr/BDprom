#!/usr/bin/env Rscript

# parse assembly metadata for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "tidyverse",
  "R.utils"
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
  make_option(
    c("-i", "--input_directory"),
    type = "character",
    default = NULL,
    help = "directory containing files from featurecounts",
    metavar = "character"
  ),
  make_option(
    c("-o", "--Gene_expression_hst2.tsv"),
    type = "character",
    default = NULL,
    help = "The output file name of the parsed featurecounts",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_directory) | is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("You must provide the name of the directory containing the featurecounts as input directory and a name for the output file.\n", call. = FALSE)
}

# Set the working directory
# this will be taken out once script is made executable
# setwd("P:/CQ_Beratug_Bildung_GmbH/Internship_2026/IOW_2026/BDP_IOW/data")

# check the input files containing files with the pattern "_counts.txt$"  and converts them into a vector

file_list_hst2 <- list.files(opt$indir, pattern = "_counts.txt$")

# append names to each element of the vector
names(file_list_hst2) <- gsub(pattern = "_counts.txt$", "", file_list_hst2)

# remove the comments rows (#) and join each row (gene name) to the 7th column (counts)
Gene_expression_hst2 <- lapply(
  names(file_list_hst2),
  function(x) {
    tmp <- read.table(file_list_hst2[x], h = T, sep = "\t", row.names = 1)[, 6, drop = F]
    colnames(tmp) <- x
    return(tmp)
  }
) |>
  bind_cols() |>
  rownames_to_column(var = "Ensembl_gene_id")

# write output
write.table(
  Gene_expression_hst2,
  file = opt$output,
  quote = T,
  sep = ",",
  row.names = F,
  col.names = T
)
