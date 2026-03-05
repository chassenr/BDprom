visualize_candidate_gene_pairs <- function(
    ge,
    col_pal = rainbow,
    out_file = "bdprom_plot_pairs.pdf"
) {
  
  study.col <- col_pal(length(ge$tpm))
  names(study.col) <- names(ge$tpm)
  
  # first, get R2 of linear fit across all studies passing thresholds in previous step to sort figures
  ge_r2 <- sapply(
    unique(ge$bdp$gene_pair),
    function(i) {
      g1 <- gsub("___.*$", "", i)
      g2 <- gsub("^.*___", "", i)

      ge1_pass <- unlist(sapply(ge$tpm[ge$bdp$study[ge$bdp$gene_pair == i]], function(x) x[, g1]))
      ge2_pass <- unlist(sapply(ge$tpm[ge$bdp$study[ge$bdp$gene_pair == i]], function(x) x[, g2]))
      fit_pass <- lm(I(ge2_pass) ~ 0 + ge1_pass)
      
      return(summary(fit_pass)$r.squared)
    }
  )
  
  # then generate figures
  pdf(out_file, height = 6, width = 7, onefile = T)
  for(i in names(ge_r2)[order(ge_r2, decreasing = T)]) {
    
    g1 <- gsub("___.*$", "", i)
    g2 <- gsub("^.*___", "", i)
    plot.lim <- c(
      min(unlist(sapply(ge$tpm, function(x) if(sum(c(g1, g2) %in% colnames(x)) == 2) min(x[, c(g1, g2)])))), 
      max(unlist(sapply(ge$tpm, function(x) if(sum(c(g1, g2) %in% colnames(x)) == 2) max(x[, c(g1, g2)]))))
    )
    
    par(mar = c(3, 3, 3, 8), xpd = F)
    plot(
      0,
      0,
      xlim = plot.lim,
      ylim = plot.lim,
      type = "n",
      xlab = g1,
      ylab = g2,
      cex.lab = 0.9,
      main = i,
      mgp = c(1.5, 0.5, 0),
      axes = F
    )
    axis(1, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    axis(2, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    for(j in names(study.col)[!names(study.col) %in% ge$bdp$study[ge$bdp$gene_pair == i]]) {
      if(g1 %in% colnames(ge$tpm[[j]]) & g2 %in% colnames(ge$tpm[[j]])) {
        points(
          ge$tpm[[j]][, g1],
          ge$tpm[[j]][, g2],
          pch = 16,
          col = "grey",
          cex = 0.5
        )
      }
    }
    for(j in ge$bdp$study[ge$bdp$gene_pair == i]) {
      points(
        ge$tpm[[j]][, g1],
        ge$tpm[[j]][, g2],
        pch = 16,
        col = study.col[j]
      )
    }
    ge1_all <- unlist(sapply(ge$tpm, function(x) if(sum(c(g1, g2) %in% colnames(x)) == 2) x[, g1]))
    ge2_all <- unlist(sapply(ge$tpm, function(x) if(sum(c(g1, g2) %in% colnames(x)) == 2) x[, g2]))
    fit_all <- lm(I(ge2_all) ~ 0 + ge1_all)
    ge1_pass <- unlist(sapply(ge$tpm[ge$bdp$study[ge$bdp$gene_pair == i]], function(x) x[, g1]))
    ge2_pass <- unlist(sapply(ge$tpm[ge$bdp$study[ge$bdp$gene_pair == i]], function(x) x[, g2]))
    fit_pass <- lm(I(ge2_pass) ~ 0 + ge1_pass)
    abline(fit_all, col = "grey", lwd = 2, lty = 2)
    abline(fit_pass, col = "black", lwd = 2, lty = 1)
    box("plot", lwd = 0.5)
    text(
      par("usr")[1],
      par("usr")[4] - (par("usr")[4] - par("usr")[3]) * 0.03,
      labels = paste0("slope = ", round(fit_pass$coefficients, 2), "; R2 = ", round(summary(fit_pass)$r.squared, 2)),
      pos = 4,
      col = "black",
      cex = 0.7
    )
    text(
      par("usr")[1],
      par("usr")[4] - (par("usr")[4] - par("usr")[3]) * 0.06,
      labels = paste0("slope = ", round(fit_all$coefficients, 2), "; R2 = ", round(summary(fit_all)$r.squared, 2)),
      pos = 4,
      col = "grey",
      cex = 0.7
    )
    par(xpd = NA)
    legend(
      x = par("usr")[2],
      y = par("usr")[4],
      legend = paste0(ge$bdp$study[ge$bdp$gene_pair == i], ": ", round(ge$bdp$cor[ge$bdp$gene_pair == i], 2)),
      col = study.col[ge$bdp$study[ge$bdp$gene_pair == i]],
      cex = 0.7,
      pt.cex = 1,
      pch = 16,
      bty = "n"
    )
  }
  dev.off()
}