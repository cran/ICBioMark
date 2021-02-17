#' Gene Lengths from the Ensembl Database
#'
#' @description Pre-imported length data from the Ensembl database for all genes on chromosomes 1-22, X and Y.
#'
#' @format A dataframe with three columns:
#' \describe{
#'   \item{Hugo_Symbol}{The names of all nuclear genes in humans for which ensembl entries with coding sequence lengths exist.}
#'   \item{max_cds}{The maximum coding sequence for each gene as given by the ensembl database.}
#'   \item{Chromosome}{The chromosome where each gene is located.}
#' }
#' @source \url{https://www.ensembl.org/}
"ensembl_gene_lengths"
