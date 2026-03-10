# todo:
# include progress messages

identify_candidate_gene_pairs <- function(
    gtf = list.files(path = ".", pattern = "*.gtf")[1],
    dist_max = 600,
    out_file = "bdprom_candidate_pairs.txt"
) {
  gene_coord <- read.table(gtf, h = F, sep = "\t") %>% 
    # only consider coordinates for transcript (alternatives: gene, transcript)
    filter(V3 == "exon") %>% 
    # extract gene names from gtf column 9
    mutate(id = gsub("gene_id ", "", gsub(";.*$", "", V9))) %>% 
    # group data set by gene name
    group_by(id) %>% 
    # in case a gene is represented by multiple lines in the gtf file (e.g. when selecting exon earlier)
    summarise(
      # get total gene length (important for quantification later)
      length = sum(V5 - V4 + 1),
      # get full extent of gene (min...
      min = min(V4),
      # ... and max)
      max = max(V5),
      # summarize orientation
      orientation = paste0(unique(V7), collapse = ";"),
      # get start position depending on orientation
      start = ifelse(orientation == "-", max, min),
      # get chromosome name to make sure that pairs aren't selected across multiple chromosomes
      chromosome = unique(V1)
    ) %>% 
    # order data by chromosome and smaller gene coordinate
    arrange(chromosome, min) %>% 
    # convert gene names to upper case to match diatomicbase output
    mutate(id = toupper(id)) %>% 
    # set gene names as rownames
    column_to_rownames("id") %>%
    # convert to data frame
    as.data.frame()
  
  # identify neighboring gene (the one upstream of the start position in reading direction)
  gene_coord$neighbor <- sapply(
    # for each gene name, i.e. each row in the data frame
    1:nrow(gene_coord), 
    function(x) {
      ifelse(
        # if gene at position x is in reverse orientation and gene at position x + 1 located on the same chromosome
        gene_coord$orientation[x] == "-" & gene_coord$chromosome[x + 1] == gene_coord$chromosome[x], 
        # the neighbor of gene x is gene x + 1
        rownames(gene_coord)[x + 1],
        ifelse(
          # if gene at position x is in forward orientation and gene at position x - 1 located on the same chromosome
          gene_coord$orientation[x] == "+" & gene_coord$chromosome[x - 1] == gene_coord$chromosome[x], 
          # the neighbor of gene x is gene x - 1
          rownames(gene_coord)[x - 1],
          # otherwise, there is no neighbor (start and end of chromosomes)
          NA
        )
      )
    }
  )
  
  # check for distance between neighbors
  candidate_list <- lapply(
    # for each gene name, i.e. each row in the data frame
    1:nrow(gene_coord),
    function(x) {
      if(
        # there is a neighbor for gene x
        !is.na(gene_coord$neighbor[x]) &
        # if gene x and its neighbor are in opposite orientation
        gene_coord$orientation[x] != gene_coord[gene_coord$neighbor[x], "orientation"] & 
        # if the distance between the neighbors is without the limits defined earlier
        abs(gene_coord$start[x] - gene_coord[gene_coord$neighbor[x], "start"]) <= dist_max
      ) {
        # keep the gene (and its neighbor) as candidates
        gene_coord[x, ]
      }
    }
  )
  # remove empty elements (i.e. non candidate gene)
  candidate_list <- candidate_list[!sapply(candidate_list, is.null)]
  candidate_pairs <- do.call("rbind", candidate_list)
  print(paste0("There are ", nrow(candidate_pairs) / 2, " candidate gene pairs."))
  
  # write cadidate genes to file to file, but also return R object
  write.table(candidate_pairs, out_file, quote = F, sep = "\t")
  out <- list(gene_coord, candidate_pairs)
  names(out) <- c("gene_coord", "gene_pairs")
  return(out)
}