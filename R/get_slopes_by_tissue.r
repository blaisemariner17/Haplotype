#'Get the reads from the recessive allele for all genes in an individual for all samples
#'
#' @param index index used to run in parallel-- the element of the list in the list_gtex_max/list_gtex_min that is to be used
#' @param list_gtex_max max reads allele by gene
#' @param list_gtex_min min reads allele by gene
#' @return Function returns matrix.
#' @export get_slopes_by_tissue

get_slopes_by_tissue <- function(index, list_gtex_max = list_gtex_max, list_gtex_min = list_gtex_min ) {
  GTEX_counts_reduced_max <- list_gtex_max[[index]]
  GTEX_counts_reduced_min <- list_gtex_min[[index]]
  result <- data.frame("ENS" = c(), "HGNC" = c(), "Tissue" = c(), "Lin_Reg_slope" = c())

  for (ens_ in rownames(GTEX_counts_reduced_max)) {
    ens <- gsub("\\..*","",ens_)
    biomart_inq <- biomaRt::getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                            "hgnc_symbol", "description"),values=ens,mart= mart)

    if (is.null(nrow(biomart_inq)) == T) {biomart_inq <- ens} else if (length(biomart_inq) == 0) {biomart_inq <- ens} else {
      hgnc <- biomart_inq$hgnc_symbol
    }

    if (length(hgnc) == 0) {hgnc <- ens}
    if (length(hgnc) > 1) {hgnc <- hgnc[1]}
    if (is.na(hgnc) == T) {hgnc <- ens}

    reads_oi_max <- GTEX_counts_reduced_max[grepl(ens, rownames(GTEX_counts_reduced_max)),]
    reads_oi_min <- GTEX_counts_reduced_min[grepl(ens, rownames(GTEX_counts_reduced_min)),]

    plotting <- data.frame(cbind(t(reads_oi_max), t(reads_oi_min)))

    if (ncol(na.omit(plotting)) > 0){
      colnames(plotting) <- c("reads_allele_max", "reads_allele_min")
      for (tissue in unique(metaData$tissue)) {
        samples_oi <- rownames(metaData)[metaData$tissue == tissue]
        metaData_oi <- metaData[metaData$tissue == tissue, ]

        plotting_ <- plotting[rownames(plotting) %in% samples_oi,]
        if (nrow(na.omit(plotting_)) > 0){
          n = nrow(plotting_)
          lm_ <- lm(reads_allele_min~reads_allele_max,data=plotting_)
          cef <- coef(lm_)
          slope <- as.numeric(cef["reads_allele_max"])
          result_ <- data.frame("ENS" = ens, "HGNC" = hgnc, "Tissue" = tissue, "N" = n, "Lin_Reg_slope" = slope)
          result <- rbind(result, result_)
        }
      }
    }
  }
  return(result)
}
