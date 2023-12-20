#### This is a script that is used to generate the data and interpretations from GTEx haplotype matrix in Mariner et al 2023
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
  "ggrepel"
)

cores_ <- 4

lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),plot.title = element_text(family = "sans", size = 12, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"), axis.text = element_text(family = "sans", size = 18, color="black"),axis.title.y = element_markdown(family = "sans", size = 20), axis.title.x = element_markdown(family = "sans", size = 20),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

# read in the data and metadata
metaData <- read.csv("../GTEx_data_metaData/metaData.csv", row.names = 1)
metaData$Age_num <- as.numeric(substr(metaData$Age,1,nchar(metaData$Age)-3))
GTEX_allelecounts <- read_rds("../GTEx_data_metaData/GTEX_countdate_byallele.rds")
if ("X" %in% colnames(GTEX_allelecounts)){
  rownames(GTEX_allelecounts) <- GTEX_allelecounts$X
  GTEX_allelecounts <- GTEX_allelecounts[,-1]
}

#biomaRt
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# seperate the gtex haplotype data into two dataframes from A reads and B reads
GTEX_allelecounts_A <- GTEX_allelecounts[grepl(".A", rownames(GTEX_allelecounts)),]
GTEX_allelecounts_B <- GTEX_allelecounts[grepl(".B", rownames(GTEX_allelecounts)),]
#rename the rownames of both of these matrices
rownames(GTEX_allelecounts_A) = substr(rownames(GTEX_allelecounts_A),1,nchar(rownames(GTEX_allelecounts_A))-2)
rownames(GTEX_allelecounts_B) = substr(rownames(GTEX_allelecounts_B),1,nchar(rownames(GTEX_allelecounts_B))-2)
rm(GTEX_allelecounts)
#generate a new matrix of total counts
GTEX_counts <- GTEX_allelecounts_A + GTEX_allelecounts_B
#time to filter genes for poor reads and samples with low number of reads
keep <- (rowSums(GTEX_counts >= 100) >= 50)
GTEX_counts_reduced <- GTEX_counts[keep,]
GTEX_allelecounts_A <- GTEX_allelecounts_A[keep,]
GTEX_allelecounts_B <- GTEX_allelecounts_B[keep,]
rm(GTEX_counts)
keep <- (colSums(GTEX_counts_reduced) >= 1000000)
GTEX_counts_reduced <- GTEX_counts_reduced[,keep]
GTEX_allelecounts_A <- GTEX_allelecounts_A[,keep]
GTEX_allelecounts_B <- GTEX_allelecounts_B[,keep]
#trim the metaData to reflect the filtering
metaData <- metaData[rownames(metaData) %in% colnames(GTEX_counts_reduced),]

#### Now we can start working with the data
library("Haplotype")



