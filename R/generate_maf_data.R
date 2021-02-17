#' Generate mutation data.
#'
#' @description A function to randomly simulate an (abridged) annotated mutation file, containing information on sample of origin, gene and mutation type, as well as a dataframe of gene lengths.
#'
#' @param n_samples (numeric)
#' The number of samples to generate mutation data for - each will have a unique value in the 'Tumor_Sample_Barcode'
#' column of the simulated MAF table. Note that if no mutations are simulated for an example, they will not appear in the table.
#' @param n_genes (numeric)
#' The number of genes to generate mutation data for - each will have a unique value in the 'Hugo_Symbol' column
#' of the simulated MAF table. A length will also be generated for each gene, and stored in the table 'gene_lengths'.
#' @param mut_types (numeric)
#' A vector of positive values giving the relative average abundance of each mutation type. The names of each mutation
#' type are stored in the names attribute of the vector, and will form the entries of the column 'Variant_Classification'
#' in the output MAF table.
#' @param data_dist (function)
#' Directly provide the probability distribution of mutations, as a function on n_samples, n_genes, mut_types, and gene_lengths.
#' @param sample_rates (numeric)
#' Directly provide sample-specific rates.
#' @param gene_rates (numeric)
#' Directly provide gene-specific rates.
#' @param gene_lengths (numeric)
#' Directly provide gene lengths, in the form of a vector of numerics with names attribute corresponding to gene names.
#' @param sample_rates_dist (function)
#' Directly provide the distribution of sample-specific rates, as a function of the number of samples.
#' @param gene_rates_dist (function)
#' Directly provide the distribution of gene-specific rates, as a function of the number of genes.
#' @param gene_lengths_dist (function)
#' Directly provide the distribution of gene lengths, as a function of the number of genes.
#' @param bmr_genes_prop (numeric)
#' The proportion of genes that follow the background mutation rate. If specified (as is automatic), this proportion of genes
#' will have gene-specific rates equal to 1. By setting to be NULL, can avoid applying this step.
#' @param output_rates (logical)
#' If TRUE, will include the sample and gene rates in the output.
#' @param seed_id (numeric)
#' Input value for the function set.seed().
#'
#' @return
#' A list with two elements, 'maf' and 'gene_lengths'. These are (respectively):
#' * A table with three columns: 'Tumor_Sample_Barcode', 'Hugo_Symbol' and 'Variant_Classification', listing the mutations occurring in
#' the simulated example.
#' gene_lengths (dataframe)
#' * A table with two rows: 'Hugo_Symbol' and 'gene_lengths'.
#' @export
#'
#' @examples
#' # Generate some random data
#' data <- generate_maf_data(n_samples = 10, n_genes = 20)
#' # See the first rows of the maf table.
#' print(head(data$maf))
#' # See the first rows of the gene_lengths table.
#' print(head(data$gene_lengths))

generate_maf_data <- function(n_samples = 100, n_genes = 20, mut_types = NULL, data_dist = NULL,
                              sample_rates = NULL, gene_rates = NULL, gene_lengths = NULL,
                              sample_rates_dist = NULL, gene_rates_dist = NULL, gene_lengths_dist = NULL,
                              bmr_genes_prop = 0.7, output_rates = FALSE, seed_id = 1234) {

  set.seed(seed_id)

  gene_list <- paste0("GENE_", 1:n_genes)
  sample_list <- paste0("SAMPLE_", 1:n_samples)

  if(is.null(mut_types)) {
    mut_types <- c(0.6, 0.2, 0.05, 0.02,  rep(0.01, 13))
    names(mut_types) <- c('Missense_Mutation', 'Silent',
                          'Nonsense_Mutation', 'Frame_Shift_Del',
                          'Splice_Site', 'Frame_Shift_Ins',
                          'Splice_Region', '3\'Flank', '5\'Flank',
                          'Intron', 'RNA', '3\'UTR', '5\'UTR',
                          'Translation_Start_Site', 'Nonstop_Mutation',
                          'In_Frame_Ins','In_Frame_Del')
  }

  n_mut_types <- length(mut_types)

  if (is.null(gene_lengths)) {
    if (is.null(gene_lengths_dist)) {
      gene_lengths_dist <- function(n) {

        return(stats::rpois(n = n, lambda = 1000))}
    }
    gene_lengths <- gene_lengths_dist(n = n_genes)
    names(gene_lengths) <- gene_list
  }

  message("Generating data")
  if (is.null(data_dist)) {
    data_dist <- function(n_samples, n_genes, mut_types, gene_lengths, output_rates) {

      if (is.null(sample_rates)) {
        if(is.null(sample_rates_dist)) {
          sample_rates_dist <- function(n) {
            exp(stats::rnorm(n = n, mean = -8, sd = 1))
          }
        }
        sample_rates <- sample_rates_dist(n = n_samples)
        names(sample_rates) <- sample_list
      }

      if (is.null(gene_rates)) {
        if (is.null(gene_rates_dist)) {
          gene_rates_dist <- function(n) {

            return(10^stats::runif(n = n, min = -1, max = 1))}
        }
        gene_rates <- gene_rates_dist(n = n_genes)
        names(gene_rates) <- gene_list
      }

      if (!is.null(bmr_genes_prop)) {
        bmr_genes <- sample(gene_list, size = floor(bmr_genes_prop*n_genes), replace = FALSE)
        gene_rates[bmr_genes] <- 1
      }

      mutation_vector_means <- rep(sample_rates, times = n_genes * n_mut_types) *
                                rep(mut_types, each = n_samples, times = n_genes) *
                                  rep(gene_rates * gene_lengths, each = n_samples * n_mut_types)

      if (output_rates) {
        return(list(vector = stats::rpois(n = n_samples * n_genes * length(mut_types), lambda = mutation_vector_means), sample_rates = sample_rates, gene_rates = gene_rates))
      }
      else{
        return(list(vector = stats::rpois(n = n_samples * n_genes * length(mut_types), lambda = mutation_vector_means)))
      }
    }
  }

  mutation_data <- data_dist(n_samples = n_samples, n_genes = n_genes, mut_types = mut_types, gene_lengths = gene_lengths, output_rates = output_rates)
  mutation_vector <- mutation_data$vector

  message("Assembling maf")
  maf <- data.frame(Tumor_Sample_Barcode = rep("", sum(mutation_vector)),
                    Hugo_Symbol = rep("", sum(mutation_vector)),
                    Variant_Classification = rep("", sum(mutation_vector)), stringsAsFactors = FALSE)


  ### Filling out maf - could be made more efficient, but probably doesn't need to be.
  index <- 1

  pb <- utils::txtProgressBar(max = n_genes, width = 100, style = 3)
  for (gene in 1:n_genes) {
    utils::setTxtProgressBar(pb, gene)

    for (mut_type in 1:n_mut_types) {
      for (sample in 1:n_samples) {
        vector_position <- (gene - 1)*n_mut_types*n_samples + (mut_type - 1)*n_samples + sample
        mutation_count <- mutation_vector[vector_position]

        if (mutation_count > 0) {
          maf[index:(index + mutation_count - 1), "Tumor_Sample_Barcode"] <- sample_list[sample]
          maf[index:(index + mutation_count - 1), "Hugo_Symbol"] <- gene_list[gene]
          maf[index:(index + mutation_count - 1), "Variant_Classification"] <- names(mut_types)[mut_type]
          index <- index + mutation_count
        }
      }
    }
  }

  maf <- maf[sample(1:sum(mutation_vector), sum(mutation_vector), replace = FALSE),]

  gene_length_data <- data.frame(Hugo_Symbol = names(gene_lengths), max_cds = gene_lengths)

  if (output_rates) {
    return(list(maf = maf, gene_lengths = gene_length_data, rates = list(sample_rates = mutation_data$sample_rates, gene_rates = mutation_data$gene_rates, mut_types = mut_types)))
  }
  else {
    return(list(maf = maf, gene_lengths = gene_length_data))
  }
}
