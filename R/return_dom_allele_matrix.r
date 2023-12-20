#'Get the reads for the dominant allele for all genes in an individual for all samples
#'
#' @param samples_oi_A dataframe of allele A counts to be passed in
#' @param samples_oi_B dataframe of allele B counts to be passed in
#' @return Function returns matrix.
#' @export return_dom_allele_matrix

return_dom_allele_matrix <- function(samples_oi_A, samples_oi_B){
  samples_oi_max <- samples_oi_A

  for (i in 1:nrow(samples_oi_A)){

    list_ = list()
    list_[[1]] <- rowSums(samples_oi_A[i,])
    list_[[2]] <- rowSums(samples_oi_B[i,])

    oi <- grep(max(list_[[1]], list_[[2]]),
               list_)

    if (length(oi) > 1) {next}

    if (oi == 2) {
      samples_oi_max[i,] <- samples_oi_B[i,]
    }
  }
  return(samples_oi_max)
}
