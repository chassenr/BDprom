# todo:
# include progress messages

confirm_candidate_gene_pairs <- function(
    candidate_pairs,
    prior_out = NULL,
    ge_dir = "diatomicbase_data",
    cor_min = 0.5,
    slope_diff = 2,
    median_diff = 1.5,
    detect_f = 0.5,
    out_file = "bdprom_ge_pairs.txt"
) {
  
  if(is.null(prior_out)) {
    # get list of all gene expression data files
    gene_exp_list <- list.files(path = ge_dir, pattern = "*.csv")
    
    # read in gene expression data
    gene_exp_df <- lapply(
      gene_exp_list,
      function(x) {
        tmp <- read.table(paste0(ge_dir, "/", x), h = T, sep = ",") %>% 
          # remove character columns that are not further needed
          select(-User_ID.sorted.by.SD, -symbol) %>% 
          # use rownames for the gene names instead of a separate column to ensure that the data frame contains only numric data
          column_to_rownames("Ensembl_gene_id")
        # remove genes without any expression
        tmp <- tmp[rowSums(tmp) > 0, ]
      }
    )
    names(gene_exp_df) <- gsub("^.*PRJ", "PRJ", gsub(".csv$", "", gene_exp_list))
    
    # to calculate gene expression in TPM properly, gene coordinates have to be read for exons earlier
    # in the case of Phaeodactylum tricornutum this does not make a difference for the identification of candidate gene pairs
    
    # calculate reads per kilobase
    gene_exp_rpk <- lapply(
      gene_exp_df,
      function(x) {
        tmp <- candidate_pairs$gene_coord[rownames(x), "length"]
        apply(x, 2, function(y) y/(tmp / 1000))
      }
    )
    
    # calculate TPM from RPK
    gene_exp_tpm <- lapply(
      gene_exp_rpk,
      function(x) {
        data.frame(t(prop.table(x, 2) * 1000000))  
      }
    )
    
    # get gene pair combinations from list of candidates
    gene_pairs <- unique(data.frame(t(apply(data.frame(rownames(candidate_pairs$gene_pairs), candidate_pairs$gene_pairs$neighbor), 1, sort))))
    rownames(gene_pairs) <- apply(gene_pairs, 1, paste, collapse = "___")
    colnames(gene_pairs) <- c("gene_1", "gene_2")
    
    # calculate summary stats, correlation and slope for candidate pairs
    gene_exp_pairs <- map_dfr(
      names(gene_exp_tpm), 
      function(x) {
        # only look at pairs where both genes are expressed in the study 
        out <- gene_pairs
        out$pair_expressed <- apply(gene_pairs, 1, function(y) sum(colnames(gene_exp_tpm[[x]]) %in% y) == 2)
        out <- out[out$pair_expressed, ]
        out$pair_expressed <- NULL
        
        # pearson correlation coefficient
        out$cor <- sapply(1:nrow(out), function(y) cor(gene_exp_tpm[[x]][, out$gene_1[y]], gene_exp_tpm[[x]][, out$gene_2[y]]))
        # slope and R2 of linear regression (this assumed that the first gene in the pair will be x)
        fit <- lapply(1:nrow(out), function(y) lm(I(gene_exp_tpm[[x]][, out$gene_2[y]]) ~ 0 + gene_exp_tpm[[x]][, out$gene_1[y]]))
        out$slope <- sapply(fit, function(y) y$coefficients)
        out$r2 <- sapply(fit, function(y) summary(y)$r.squared)
        # add min, median and max gene expression
        out$min_1 <- apply(gene_exp_tpm[[x]], 2, min)[out$gene_1]
        out$min_2 <- apply(gene_exp_tpm[[x]], 2, min)[out$gene_2]
        out$median_1 <- apply(gene_exp_tpm[[x]], 2, median)[out$gene_1]
        out$median_2 <- apply(gene_exp_tpm[[x]], 2, median)[out$gene_2]
        out$max_1 <- apply(gene_exp_tpm[[x]], 2, max)[out$gene_1]
        out$max_2 <- apply(gene_exp_tpm[[x]], 2, max)[out$gene_2]
        # add study name to output table
        out$study <- x
        # include rowname as column in output table
        out <- out %>% rownames_to_column(var = "gene_pair")
        # return output table
        return(out)
      }
    )
  } else {
    # if the analysis has been run before, use existing gene expression data
    gene_exp_tpm <- prior_out$tpm
    gene_exp_pairs <- prior_out$ge
  }
  
  # filter output
  gene_exp_filt <- gene_exp_pairs
  if(!is.null(cor_min)) {
    gene_exp_filt <- gene_exp_filt %>% filter(cor >= cor_min)
  }
  if(!is.null(slope_diff)) {
    gene_exp_filt <- gene_exp_filt %>% filter(slope <= slope_diff & slope >= 1/slope_diff)
  }
  if(!is.null(median_diff)) {
    gene_exp_filt <- gene_exp_filt %>% filter(abs(median_1 - median_2) / ((median_1 + median_2) / 2) <= median_diff)
  }
  
  # consensus selection across all studies
  if(!is.null(detect_f)) {
    select_pairs <- names(which(table(gene_exp_filt$gene_pair) / length(gene_exp_tpm) >= detect_f))
  } else {
    select_pairs <- names(table(gene_exp_filt$gene_pair))
  }
  # subset table to only those gene pairs passing filter
  bidir_promoter_pairs <- gene_exp_filt[gene_exp_filt$gene_pair %in% select_pairs, ]
  # order by gene pair
  bidir_promoter_pairs <- bidir_promoter_pairs[order(bidir_promoter_pairs$gene_pair), ]
  
  # write output
  write.table(bidir_promoter_pairs, out_file, quote = F, sep = "\t", row.names = F)
  
  # return R list with relevant objects for visualization
  out_list <- list(gene_exp_tpm, gene_exp_pairs, bidir_promoter_pairs)
  names(out_list) <- c("tpm", "ge", "bdp")
  return(out_list)
}