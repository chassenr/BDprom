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
      # get end position depending on orientation
      end = ifelse(orientation == "-", min, max),
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
  
  # identify neighboring genes (the one upstream of the start position in reading direction)
  # already filter to pairs with opposite orientation
  gene_coord$neighbor <- sapply(
    # for each gene name, i.e. each row in the data frame
    1:nrow(gene_coord), 
    function(x) {
      # if gene at position x is in reverse orientation
      if(gene_coord$orientation[x] == "-") {
        # only consider genes on the same chromosome wich are located downstream (based on order in GTF)
        tmp <- gene_coord[(x + 1):nrow(gene_coord), ] %>% 
          filter(chromosome == gene_coord$chromosome[x])
        # select the first gene that has a start position higher than x and is in forward orientation
        suppressWarnings(return(rownames(tmp)[min(which(tmp$start > gene_coord$start[x] & tmp$orientation == "+"))]))
      }
      # if gene at position x is in forward orientation
      if(gene_coord$orientation[x] == "+") {
        # only consider genes on the same chromosome wich are located upstream (based on order in GTF)
        tmp <- gene_coord[1:(x - 1), ] %>% 
          filter(chromosome == gene_coord$chromosome[x])
        # select the last gene that has a start position lower than x and is in reverse orientation
        suppressWarnings(return(rownames(tmp)[max(which(tmp$start < gene_coord$start[x] & tmp$orientation == "-"))]))
      }
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
        # gene pairs are bidirectional (no need to check if neighbor of neighbor is NA - that cannot happen)
        gene_coord[gene_coord$neighbor[x], "neighbor"] == rownames(gene_coord)[x] &
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