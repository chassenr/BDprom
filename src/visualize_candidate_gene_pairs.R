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
  pdf(out_file, height = 6, width = 8, onefile = T)
  for(i in names(ge_r2)[order(ge_r2, decreasing = T)]) {
    
    g1 <- gsub("___.*$", "", i)
    g2 <- gsub("^.*___", "", i)
    plot.lim <- c(
      min(unlist(sapply(ge$tpm, function(x) if(sum(c(g1, g2) %in% colnames(x)) == 2) min(x[, c(g1, g2)])))), 
      max(unlist(sapply(ge$tpm, function(x) if(sum(c(g1, g2) %in% colnames(x)) == 2) max(x[, c(g1, g2)]))))
    )
    
    layout(mat = matrix(c(1, 1, 1, 2, 3, 4), nrow = 3, ncol = 2, byrow = F), widths = c(1.5, 1))
    par(mar = c(3, 3, 3, 1), yaxs = "r", xpd = F, lwd = 0.1)
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
          col = "grey80",
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
    abline(fit_all, col = "grey80", lwd = 2, lty = 2)
    abline(fit_pass, col = "grey20", lwd = 2, lty = 1)
    box("plot", lwd = 0.5)
    text(
      par("usr")[1],
      par("usr")[4] - (par("usr")[4] - par("usr")[3]) * 0.03,
      labels = paste0("slope = ", round(fit_pass$coefficients, 2), "; R2 = ", round(summary(fit_pass)$r.squared, 2)),
      pos = 4,
      col = "grey20",
      cex = 0.7
    )
    text(
      par("usr")[1],
      par("usr")[4] - (par("usr")[4] - par("usr")[3]) * 0.06,
      labels = paste0("slope = ", round(fit_all$coefficients, 2), "; R2 = ", round(summary(fit_all)$r.squared, 2)),
      pos = 4,
      col = "grey80",
      cex = 0.7
    )
    
    par(mar = c(3, 2, 3, 1), yaxs = "i", xpd = NA)
    # correlation per study
    hist(
      ge$ge$cor[ge$ge$gene_pair == i],
      breaks = seq(-1, 1, 0.1),
      col = "grey80",
      axes = F,
      main = "Correlation coefficient",
      ylab = "", 
      xlab = ""
    )
    axis(1, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    axis(2, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    hist(
      ge$bdp$cor[ge$bdp$gene_pair == i],
      breaks = seq(-1, 1, 0.1),
      col = "grey20",
      add = T
    )
    text(
      par("usr")[1],
      par("usr")[4] * 1.05,
      labels = paste0(sum(ge$bdp$gene_pair == i), " out of ", sum(ge$ge$gene_pair == i), " studies passing thresholds."),
      pos = 4,
      col = "grey20",
      cex = 0.7
    )
    # R2 per study
    hist(
      ge$ge$r2[ge$ge$gene_pair == i],
      breaks = seq(0, 1, 0.05),
      col = "grey80",
      axes = F,
      main = "R squared",
      ylab = "", 
      xlab = ""
    )
    axis(1, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    axis(2, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    hist(
      ge$bdp$r2[ge$bdp$gene_pair == i],
      breaks = seq(0, 1, 0.05),
      col = "grey20",
      add = T,
      lwd = 0.1
    )
    # slope per study
    hist(
      ge$ge$slope[ge$ge$gene_pair == i],
      breaks = seq(min(ge$ge$slope[ge$ge$gene_pair == i]), max(ge$ge$slope[ge$ge$gene_pair == i]), length.out = 20),
      col = "grey80",
      axes = F,
      main = "Slope",
      ylab = "", 
      xlab = ""
    )
    axis(1, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    axis(2, lwd = 0.1, cex.axis = 0.7, mgp = c(1.5, 0.5, 0))
    hist(
      ge$bdp$slope[ge$bdp$gene_pair == i],
      breaks = seq(min(ge$ge$slope[ge$ge$gene_pair == i]), max(ge$ge$slope[ge$ge$gene_pair == i]), length.out = 20),
      col = "grey20",
      add = T
    )
  }
  dev.off()
}