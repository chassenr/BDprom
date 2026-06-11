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

# Some functions for message output
msg <- function(X) {
  cat(crayon::white(paste0("[", format(Sys.time(), "%T"), "]")), X)
}


#####################################################################
### Reading command line options ####

# define command line options
option_list <- list(
  make_option(
    c("-e", "--ena_search_out"),
    type = "character",
    default = NULL,
    help = "import TSV search output obtained from ENA portal API",
    metavar = "character"
  ),
  make_option(
    c("-c", "--min_readcounts"),
    type = "integer",
    default = 5000000,
    help = "set a minimum read count value",
    metavar = "integer"
  ),
  make_option(
    c("-r", "--min_runs"),
    type = "integer",
    default = 5,
    help = "set minimum number of runs",
    metavar = "integer"
  ),
  make_option(
    c("-s", "--seq_tech"),
    type = "character",
    default = "ILLUMINA,BGISEQ,DNBSEQ",
    help = "set the requires sequencing technology",
    metavar = "character"
  ),
  make_option(
    c("-p", "--pubmed"),
    type = "logical",
    default = FALSE,
    help = "retain pubmed entries only",
    metavar = "logical"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "output file name",
    metavar = "character"
  )
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$ena_search_out) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("You must provide the name of the ENA search output file in tsv format as input and a name for the output file.\n", call. = FALSE)
}


# Set the working directory
# this will be taken out once script is made executable
# setwd("P:/CQ_Beratug_Bildung_GmbH/Internship_2026/IOW_2026/BDP_IOW/data")

#  parse vector of sequencing platform names from opt
seq_tech <- strsplit(opt$seq_tech, ",")[[1]]

# reading data
ENA_curl <- read.table(file = opt$ena_search_out, comment.char = "", quote = "", sep = "\t", h = T) |>
  filter(
    run_accession != "", 
    experiment_accession != "",
    sample_accession != "",
    study_accession != ""
  )
# head(ENA_curl)
# colnames(ENA_curl)

# Getting Xrefs for the study accessions
# Any unique study accessions?

msg(
  paste0(
    "There are ",
    nrow(ENA_curl),
    " runs from ",
    length(unique(ENA_curl$study_accession)),
    " ENA studies. The total data volume is ",
    hsize(sum(sapply(strsplit(ENA_curl$fastq_bytes, ";"), function(x) sum(as.numeric(x))))),
    ".\n"
  )
)


# Xref links retrieval
# test <- read.table("https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession=PRJEB11970", header =TRUE, sep = "\t")

combined_xrefs <- lapply(
  unique(ENA_curl$study_accession),
  function(study_acc) {
    read.table(
      url(paste0("https://www.ebi.ac.uk/ena/xref/rest/tsv/search?accession=", study_acc)),
      h = F,
      sep = "\t",
      col.names = c("DB", "DB_acc", "DB_sec", "DB_url", "study", "study_accession", "study_sec", "study_url"),
      colClasses = rep("character", 8)
    )
  }
) |>
  # combine everything into one output
  bind_rows() |>
  # Filter in only entries from EuropePMC and PubMed entries only
  filter(DB %in% c("EuropePMC", "PubMed")) |>
  # Remove redundancies. Create separate columns to separate EuropePMC from PubMed
  mutate(study_primary = ifelse(grepl("^P", study_accession), study_accession, study_sec)) |>
  select(DB, DB_acc, study_primary) |>
  pivot_wider(
    names_from = DB,
    id_cols = study_primary,
    values_from = DB_acc,
    values_fn = function(x) paste(unique(x), collapse = ";")
  )

msg(paste0(nrow(combined_xrefs), " of those ENA studies have associated scientific articles.\n"))


# Now, combine the combined_xref_PubMed table to the very original ENA_curl table that we had obtained from the ENA database. #left_join()

ENA_xrefs <- left_join(
  ENA_curl,
  combined_xrefs,
  by = c("study_accession" = "study_primary")
)

# colnames(ENA_xrefs_filtered)
# nrow(ENA_xrefs)

#### filter metadata
ENA_xrefs_filtered <- ENA_xrefs |>
  # remove runs without correct data i.e. filter out entries that lack downloadable .fastq files
  filter(
    !is.na(fastq_ftp) & fastq_ftp != "",
    # remove shallowly sequenced runs
    read_count >= opt$min_readcounts
  ) |>
  # remove fewer than 5 replicates/runs i.e. >= 5, first group by study accession (by rows), then filter and count the runs that meet the specified threshold.

  # by the type of the instrument platform
  # since the allowed platforms are listed in a vector, use %in% to go inside the vector

  filter(instrument_platform %in% seq_tech) |>
  # for the PE reads, it can happen that the entries just indicates PE (due to poor documentation and reportin.     #For this, we need to split the fastq_ftp files that have entries separated by semicolon (the are the 2 files of the PAIRED ENDS).
  # create a new column to save these the counts (from the function length) of the splitted .fastq files.
  mutate(n_files = sapply(strsplit(fastq_ftp, ";"), length)) |>
  # filter only those that are recorded as paired and have 2 .fast files from the spli, AND those recorded as SINGLE and have 1 .fastq file files
  filter((library_layout == "PAIRED" & n_files == 2) | (library_layout == "SINGLE" & n_files == 1)) |>
  # count number of runs per study and filter by minimum number of runs
  group_by(study_accession) |>
  filter(n() >= opt$min_runs) |>
  ungroup()


msg(
  paste0(
    "After filtering there are now ",
    nrow(ENA_xrefs_filtered),
    " runs from ",
    length(unique(ENA_xrefs_filtered$study_accession)),
    " ENA studies. The total data volume is ",
    hsize(sum(sapply(strsplit(ENA_xrefs_filtered$fastq_bytes, ";"), function(x) sum(as.numeric(x))))),
    ". Reads counts range from ",
    min(ENA_xrefs_filtered$read_count),
    " to ",
    max(ENA_xrefs_filtered$read_count), ".\n"
  )
)

# Retain only entries whose Pumbmed articles are available
if (opt$pubmed) {
  ENA_xrefs_filtered <- ENA_xrefs_filtered |> filter(!is.na(PubMed) | !is.na(EuropePMC))

  msg(
    paste0(
      "Removing studies without an associated scientific article resulted in ",
      nrow(ENA_xrefs_filtered),
      " runs from ",
      length(unique(ENA_xrefs_filtered$study_accession)),
      " ENA studies being retained. The total data volume is ",
      hsize(sum(sapply(strsplit(ENA_xrefs_filtered$fastq_bytes, ";"), function(x) sum(as.numeric(x))))),
      ". Reads counts range from ",
      min(ENA_xrefs_filtered$read_count),
      " to ",
      max(ENA_xrefs_filtered$read_count), ".\n"
    )
  )
}


## format download links

#ENA_xrefs_downloadlinks <- ENA_xrefs_filtered |>
#  separate(
#    col = fastq_ftp,
#    into = c("r1_download", "r2_download"),
#    sep = ";",
#    remove = F,
#    convert = F,
#    fill = "right"
#  )

# write output
write.table(
  ENA_xrefs_filtered,
  opt$output,
  quote = F,
  sep = "\t",
  row.names = F
)
