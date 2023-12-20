#'Get the reads for the recessive allele for all genes in an individual for all samples
#'
#' @param ind individuals of interest
#' @param GTEX_allelecounts_A allelecounts A
#' @param GTEX_allelecounts_B allelecounts B
#' @param metaData the metadata which tells us which individuals belong to which columns in the allele data
#' @return Function returns matrix of the recessive allele for each individual
#' @export get_min_by_sample

get_min_by_sample <- function(ind, metaData = metaData, GTEX_allelecounts_A= GTEX_allelecounts_A, GTEX_allelecounts_B = GTEX_allelecounts_B) {
  samples_oi <- rownames(metaData)[metaData$individual == ind]
  samples_oi_A <- GTEX_allelecounts_A[,colnames(GTEX_allelecounts_A) %in% samples_oi]
  samples_oi_B <- GTEX_allelecounts_B[,colnames(GTEX_allelecounts_B) %in% samples_oi]
  samples_oi_max <- return_rec_allele_matrix(samples_oi_A, samples_oi_B)
  return(samples_oi_max)
}
