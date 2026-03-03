library(ggplot2)

file.names <- c("age.on.ldl", "sex.on.testosterone", "sex.on.urate", 
                "statin.on.ldl", "statin.on.a1c", "statin.on.glucose", 
                "smoke_status.on.bmi")
env.names <- sapply(strsplit(file.names, split = "\\."), function(x) x[1])
tissues <- unlist(data.table::fread("data/tissues.txt", header = FALSE))
names(tissues) <- NULL
res <- list()

for (i in seq_along(file.names)) {
  res.tmp <- list()
  file.name <- file.names[i]
  env.name <- env.names[i]
  
  for (tissue in tissues) {
    dat <- data.table::fread(
      sprintf("data/stats/ukb.txewas.%s.Tiss_spec.%s.glmsve", file.name, 
              tissue),
      sep = "\t", quote = FALSE, header = TRUE, data.table = FALSE)
    blacklist <- unlist(data.table::fread(
      sprintf("data/blacklist/gene_blacklist.pval_0.05.Tiss_spec.%s.txt", 
              tissue),
      header = FALSE))
    
    dat <- cbind(
      dat[dat$TYPE == paste("ADD", env.name, sep = "x"), c(1, 5, 6)],
      dat[dat$TYPE == "ADD", 5:6])
    colnames(dat) <- c("ID", "Z_GxE", "P_GxE", "Z_G", "P_G")
    dat <- dat[!(dat$ID %in% blacklist),]
    res.tmp[[length(res.tmp) + 1]] <- dat
  }
  
  res.tmp <- do.call("rbind", res.tmp)
  res.tmp <- res.tmp[order(abs(res.tmp$Z_GxE), decreasing = TRUE),]
  res.tmp <- res.tmp[!duplicated(res.tmp$ID),]
  
  signif.genes <- data.table::fread(
    sprintf("data/egenes/eGenes.glmsve.%s.Tiss_spec.tsv", file.name),
    sep = "\t", header = TRUE, data.table = FALSE)
  signif.genes <- signif.genes[, 1]
  signif.genes <- signif.genes[signif.genes != "testGene"]
  
  res.tmp <- res.tmp[res.tmp$ID %in% signif.genes,]
  res.tmp <- res.tmp[res.tmp$P_G < 0.05,]
  
  sgn.prod <- sign(res.tmp$Z_GxE * res.tmp$Z_G)
  npos <- sum(sgn.prod == 1)
  nneg <- sum(sgn.prod == -1)
  
  if (npos > nneg) {
    res[[i]] <- npos/length(sgn.prod)
  } else {
    res[[i]] <- nneg/length(sgn.prod)
  }
}

res <- unlist(res)
dat <- data.frame(name = c("Age on LDL", "Sex on testosterone", "Sex on urate", 
                           "Statins on LDL", "Statins on A1c", 
                           "Statins on glucose", "Smoking status on BMI"),  
                  value = res)
dat <- dat[order(dat$value, decreasing = TRUE),]
vals <- dat$value
lbls <- dat$name
dat <- do.call(
  "rbind", 
  apply(dat, 1, 
        function(x) data.frame(name = rep(x[1], 2), 
                               condition = c("plus", "minus"), 
                               value = c(as.numeric(x[2]), 1-as.numeric(x[2])))))
dat$name <- factor(dat$name, levels = lbls)

pdf("plot_frac_sign_consistent.pdf")
ggplot(dat, aes(x = name, y = value, fill = condition)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("royalblue", 'firebrick2'), guide = "none") +
  coord_flip() +
  geom_text(aes(x = name, 
                y = value, 
                label = c(rbind(sprintf("%.2f", vals), rep("", length(lbls)))), 
                group = condition),
            color = "white", size = 5.5, hjust = 1.3) +
  ylab(expression(paste("Fraction of sign-consistent ", italic("G x E")))) +
  theme_set(theme_bw(base_size = 20)) + 
  theme(axis.text = element_text(colour='black'),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())
dev.off()
