library(ggplot2)

name <- "sex.on.testosterone" # Set this to one of the following:
# "sex.on.testosterone", 
# "sex.on.urate", 
# "smoke_status.on.bmi", 
# "age.on.ldl", 
# "statin.on.a1c", 
# "statin.on.glucose", 
# "statin.on.ldl"

name.list <- strsplit(name, split = "\\.")[[1]]
env.name <- name.list[1]
outcome.name <- name.list[3]

tissues <- unlist(data.table::fread("data/tissues.txt", header = FALSE))
names(tissues) <- NULL
res <- list()

for (tissue in tissues) {
  dat <- data.table::fread(
    sprintf("data/stats/ukb.txewas.%s.Tiss_spec.%s.glmsve", name, tissue),
    sep = "\t", quote = FALSE, header = TRUE, data.table = FALSE)
  blacklist <- unlist(data.table::fread(
    sprintf("data/blacklist/gene_blacklist.pval_0.05.Tiss_spec.%s.txt", 
            tissue),
    header = FALSE))
  
  dat <- cbind(dat[dat$TYPE == paste("ADD", env.name, sep = "x"), c(1, 3)], 
               dat[dat$TYPE == "ADD", 3],
               dat[dat$TYPE == env.name, 3])
  colnames(dat) <- c("ID", "Z_GxE", "Z_G", "Z_E")
  dat <- dat[!(dat$ID %in% blacklist),]
  res[[length(res) + 1]] <- dat
}

res <- do.call("rbind", res)
res <- res[order(abs(res$Z_GxE), decreasing = TRUE),]
res <- res[!duplicated(res$ID),]

signif.genes <- data.table::fread(
  sprintf("data/egenes/eGenes.glmsve.%s.Tiss_spec.tsv", name),
  sep = "\t", header = TRUE, data.table = FALSE)
signif.genes <- signif.genes[, 1]
signif.genes <- signif.genes[signif.genes != "testGene"]

res <- res[res$ID %in% signif.genes,]
res$col <- factor(ifelse(res[, 2]*res[, 3] > 0, 1, 0), levels = c(0, 1))
betaE <- round(res[1, 4], digits = 2)

if (sum(res$col == 0) > sum(res$col == 1)) {
  cols <- c('firebrick2', "royalblue")
} else {
  cols <- c("royalblue", 'firebrick2')
}

pdf(sprintf("plot_main_vs_inter_zscores_%s_on_%s.pdf", env.name, outcome.name))
ggplot(res, aes(x = Z_G, y = Z_GxE, color = col)) + 
  geom_point() + 
  xlab(expression(beta)) +
  ylab(expression(gamma)) + 
  theme_set(theme_bw(base_size = 20)) + 
  theme(axis.text = element_text(colour='black'),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
        plot.margin = unit(c(.4, .6, .4, .1), "cm")) + 
  scale_color_manual(values = c("0" = cols[1], "1" = cols[2]), 
                     guide = "none") +
  geom_hline(yintercept = 0, color = "gray38") + 
  geom_vline(xintercept = 0, color = "gray38") + 
  annotate("text", x=-Inf, y=-Inf, hjust=-.3, vjust=-1.3, size = 5.5, 
           label = paste0("alpha == ", betaE), parse = TRUE) +
  ggtitle(paste(stringr::str_to_title(env.name), "on", outcome.name))
dev.off()
