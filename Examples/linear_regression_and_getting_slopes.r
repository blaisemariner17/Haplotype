#### This is a script that is used to generate the data and interpretations from GTEx haplotype matrix in Mariner et al 2024
# all questions should be directed towards blaisemariner17@gmail.com

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

library_list <- c(
  "tidyverse",
  "bsseq",
  "Biostrings",
  "parallel",
  "plyr",
  "foreach",
  "edgeR",
  "ggplot2",
  "ggtext",
  "ggExtra",
  "biomaRt",
  "ggrepel",
  "patchwork"
)

cores_ <- 4

lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),plot.title = element_text(family = "sans", size = 12, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"), axis.text = element_text(family = "sans", size = 18, color="black"),axis.title.y = element_markdown(family = "sans", size = 20), axis.title.x = element_markdown(family = "sans", size = 20),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

# read in the data and metadata
metaData <- read.csv("../GTEx_data_metaData/metaData.csv", row.names = 1)
metaData$Age_num <- as.numeric(substr(metaData$Age,1,nchar(metaData$Age)-3))

GTEX_counts_reduced_max <- read_rds("max_allele_matrix_by_ind.rds")
GTEX_counts_reduced_min <- read_rds("min_allele_matrix_by_ind.rds")

n <- 100
nr <- nrow(GTEX_counts_reduced_max)
list_gtex_max <- split(GTEX_counts_reduced_max, gl(ceiling(nr/n), n, nr))
list_gtex_min <- split(GTEX_counts_reduced_min, gl(ceiling(nr/n), n, nr))

INDEX <- 1:length(list_gtex_max)

time_start <- Sys.time()
print(time_start)
slopes_gene_by_tissue=parallel::mclapply(
  INDEX,
  get_slopes_by_tissue,
  list_gtex_max = list_gtex_max,
  list_gtex_min = list_gtex_min,
  mc.cores = cores_)
time_end <- Sys.time()
print(time_end - time_start)

lin_reg_res <- slopes_gene_by_tissue[[1]]
for (i in 2:length(slopes_gene_by_tissue)) {
  lin_reg_res <- rbind(lin_reg_res, slopes_gene_by_tissue[[i]])
}
lin_reg_res <- lin_reg_res[order(lin_reg_res$Lin_Reg_slope),]
head(lin_reg_res)

write_rds(x  = lin_reg_res, file = "LinReg_slopes_genes_by_tissue.rds")

plots_ <- list()
for (tissue in unique(lin_reg_res$Tissue)) {
  for_ggplot <- lin_reg_res[lin_reg_res$Tissue == tissue,]
  by_organ <- ggplot(for_ggplot, aes(x = Lin_Reg_slope)) +
    geom_histogram(fill = "black", color = "white") +
    xlab("Mean slope from linear regression") +
    ylab("Count") +
    ggtitle(paste0(tissue))+
    theme_blaise+
    # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
    scale_y_continuous(expand = c(0.01,0))+
    scale_x_continuous(expand = c(0.005,0))

  plots_[[tissue]] <- by_organ

}

pdf(file= paste0("tissue_LinRegSlope_histograms.pdf" ))
# create a 2X2 grid
par( mfrow= c(2,2) )
for (i in 1:length(plots_)) {
  (plots_[[i]])
}
dev.off()

##### but how many genes never go monoallelic vs how many are always etc.

res <- data.frame("Gene" = c(), "N_times_monoallelic" =c(), "N_times_biallelic" = c(), "N_times_total" =c())
for (gene in unique(lin_reg_res$HGNC)) {
  n_times_total <- nrow(lin_reg_res[lin_reg_res$HGNC == gene,])
  n_times_mono <- nrow(lin_reg_res[lin_reg_res$HGNC == gene & lin_reg_res$Lin_Reg_slope < 0.2,])
  n_times_bi <- nrow(lin_reg_res[lin_reg_res$HGNC == gene & lin_reg_res$Lin_Reg_slope > 0.8,])
  res_ <- data.frame("Gene" = gene, "N_times_monoallelic" = n_times_mono, "N_times_biallelic" = n_times_bi, "N_times_total" = n_times_total)

  res <- rbind(res, res_)
}

res$perc_mono <- res$N_times_monoallelic / res$N_times_total
res$perc_bi <- res$N_times_biallelic / res$N_times_total

plots_ <- list()
plot1 <- ggplot(res, aes(x = perc_mono)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("Fraction times monoallelic") +
  ylab("Count") +
  ggtitle("")+
  theme_blaise+
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite("frac_genes_monoallelic.svg", fix_text_size = F)
plot1
dev.off()

plot2 <- ggplot(res, aes(x = perc_bi)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("Fraction times biallelic") +
  ylab("Count") +
  ggtitle("")+
  theme_blaise+
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))
svglite("frac_genes_biallelic.svg", fix_text_size = F)
plot2
dev.off()
pdf(file= paste0("LinRegSlope_frac_genes.pdf" ))
# create a 2X2 grid
par( mfrow= c(2,2) )
plot1
plot2
dev.off()
