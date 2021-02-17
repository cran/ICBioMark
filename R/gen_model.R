#' Fit Generative Model
#'
#' @description A function to fit a generative model to a mutation dataset. At its heart, requires a gene_lengths dataframe (for examples of the correct format for this see the pre-loaded datasets example_maf_data$gene_lengths and ensembl_gene_lengths), and a mutation dataset. This is best supplied through the 'table' argument, and constructed via the function get_mutation_tables().
#'
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param matrix (Matrix::sparseMatrix)
#' A mutation matrix, such as produced by the function get_table_from_maf().
#' @param sample_list (character)
#' The set of samples to be modelled.
#' @param gene_list (character)
#' The set of genes to be modelled.
#' @param mut_types_list (character)
#' The set of mutation types to be modelled.
#' @param col_names (character)
#' The column names of the 'matrix' parameter.
#' @param table (list)
#' Optional parameter combining matrix, sample_list, gene_list, mut_types_list, col_names, as is produced by the function get_tables().
#' @param nlambda (numeric)
#' The length of the vector of penalty weights, passed to the function glmnet::glmnet().
#' @param n_folds (numeric)
#' The number of cross-validation folds to employ.
#' @param maxit (numeric)
#' Technical parameter passed to the function glmnet::glmnet().
#' @param seed_id (numeric)
#' Input value for the function set.seed().
#' @param progress (logical)
#' Show progress bars and text.
#' @param alt_model_type (character)
#' Used to call an alternative generative model type such as "US" (no sample-dependent parameters) or
#' "UI" (no gene/variant-type interactions).
#'
#' @return A list comprising three objects:
#' * An object 'fit', a fitted glmnet model.
#' * A table 'dev', giving average deviances for each regularisation penalty factor and cross-validation fold.
#' * An integer 's_min', the index of the regularsisation penalty minimising cross-validation deviance.
#' * A list 'names', containing the sample, gene, and mutation type information of the training data.
#' @export
#'
#' @examples
#' example_gen_model <- fit_gen_model(example_maf_data$gene_lengths, table = example_tables$train)
#' print(names(example_gen_model))


fit_gen_model <- function(gene_lengths, matrix = NULL, sample_list = NULL, gene_list = NULL, mut_types_list = NULL, col_names = NULL,
                          table = NULL, nlambda = 100, n_folds = 10, maxit = 1e9, seed_id = 1234, progress = FALSE, alt_model_type = NULL) {

  if (!is.null(alt_model_type)) {
    if (alt_model_type == "US") {
      return(fit_gen_model_unisamp(gene_lengths, matrix = NULL, sample_list = NULL, gene_list = NULL, mut_types_list = NULL, col_names = NULL,
                                   table = NULL, nlambda = 100, n_folds = 10, maxit = 1e9, seed_id = 1234, progress = FALSE))
    }

    if (alt_model_type == "UI") {
      return(fit_gen_model_uninteract(gene_lengths, matrix = NULL, sample_list = NULL, gene_list = NULL, mut_types_list = NULL, col_names = NULL,
                                      table = NULL, nlambda = 100, n_folds = 10, maxit = 1e9, seed_id = 1234, progress = FALSE))
    }
  }

  set.seed(seed_id)
  if (progress) {trace.it = 1}
  else {trace.it = 0}

  if (is.null(table) & any(is.null(matrix), is.null(sample_list), is.null(gene_list), is.null(mut_types_list), is.null(col_names))) {
    stop("If not providing a full tables object, must provide the inputs sample_list, gene_list, mut_types_list and col_names")
  }

  if(!is.null(table)) {
    matrix <- table$matrix
    sample_list <- table$sample_list
    gene_list <- table$gene_list
    mut_types_list <- table$mut_types_list
    col_names <- table$col_names
  }

  n_samples <- length(sample_list)
  n_genes <- length(gene_list)
  n_mut_types <- length(mut_types_list)

  if (any(dim(matrix) != c(n_samples, n_genes * n_mut_types))) {
    stop(paste0("Matrix has dimension ", dim(matrix), ", should have dimension ", c(n_samples, n_genes * n_mut_types)))
  }

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

  # Making output vector
  mutation_vector <- matrix
  dim(mutation_vector) <- c(n_samples * n_genes * n_mut_types, 1)

  # Constructing design matrix
  i <- 1:(n_samples*n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types), each = n_samples)
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_samples*n_mut_types)

  design_matrix <- Matrix::sparseMatrix(c(i, i), c(j, j_))

  t_i_getter <- Matrix::sparseMatrix(i = rep(1:n_samples, n_genes*n_mut_types), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_samples, n_samples*n_genes*n_mut_types))
  t_i <- as.vector(t_i_getter %*% mutation_vector)
  t_i_getter <- NULL

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  # We use weights (rather than offsets) to fit the model, as it makes glmnet more efficient.
  weights <- rep(t_s, each = n_samples, times = n_genes)* rep(t_i, times = n_mut_types*n_genes)
  weights <- weights*rep(gene_lengths[gene_list,]$max_cds, each = n_samples*n_mut_types)
  weights <- weights/mean(weights)

  weighted_observations <- ifelse(weights == 0, 0, as.vector(mutation_vector)/weights)

  # First (full dataset) run
  if (progress) {print("First glmnet run (full dataset)")}
  fit <- glmnet::glmnet(x = design_matrix, y = weighted_observations, nlambda = nlambda, weights = weights, family = "poisson", trace.it = trace.it, maxit = maxit)

  # Cross-validation
  if (progress) {print("Cross-validation:")}
  partitions <- sample(rep(1:n_folds, n_samples * n_genes * n_mut_types / n_folds), n_samples * n_genes * n_mut_types)
  dev <- matrix(0, n_folds, length(fit$lambda))

  for (fold in 1:n_folds) {
    part.design <- design_matrix[partitions != fold, ]
    part.response <- weighted_observations[partitions != fold]
    part.weights <- weights[partitions != fold]

    if (progress) {writeLines(paste("\nFitting glmnet on fold", fold))}
    part.fit <- glmnet::glmnet(part.design, y = part.response, family = "poisson", lambda = fit$lambda, weights = part.weights, trace.it = trace.it, maxit = maxit)
    part.weights <- NULL; part.design <- NULL; part.response <- NULL;

    if (progress) {print("Computing statistics")}
    ### This code should be made much simpler - it was originally written to squeeze the maximum amount out of a laptop without crashing.
    if (progress) {pb <- utils::txtProgressBar(max = 15, width = 100, style = 3)}

    part.test.design <- design_matrix[partitions == fold,]
    if (progress) {utils::setTxtProgressBar(pb, 1)}

    part.test.product <- part.test.design %*% part.fit$beta
    if (progress) {utils::setTxtProgressBar(pb, 2)}

    part.test.design <- NULL;

    part.test.product <- part.test.product + part.fit$a0
    if (progress) {utils::setTxtProgressBar(pb, 3)}

    part.fit <- NULL;

    part.test.product <- exp(part.test.product)
    if (progress) {utils::setTxtProgressBar(pb, 4)}

    part.test.weights <- Matrix::Matrix(rep(weights[partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda))
    if (progress) {utils::setTxtProgressBar(pb, 5)}

    part.test.predict <- part.test.weights * part.test.product
    if (progress) {utils::setTxtProgressBar(pb, 6)}

    part.test.product <- NULL;  part.test.weights <- NULL

    part.test.response <- Matrix::Matrix(rep(mutation_vector[1:(n_samples*n_mut_types*n_genes),][partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    if (progress) {utils::setTxtProgressBar(pb, 7)}

    part.test.residuals <- part.test.response - part.test.predict
    if (progress) {utils::setTxtProgressBar(pb, 8)}

    part.test.logpredict <- ifelse(part.test.predict == 0, 0, log(part.test.predict))
    if (progress) {utils::setTxtProgressBar(pb, 9)}

    part.test.predict <- NULL

    log_response <- mutation_vector[partitions == fold, ]
    log_response <- ifelse(log_response == 0, 0, log(log_response))
    if (progress) {utils::setTxtProgressBar(pb, 10)}

    part.test.logresponse <- Matrix::Matrix(rep(log_response, length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    if (progress) {utils::setTxtProgressBar(pb, 11)}

    log_response <- NULL

    log_div <- part.test.logresponse - part.test.logpredict
    if (progress) {utils::setTxtProgressBar(pb, 12)}

    part.test.logresponse <- NULL; part.test.logpredict <- NULL;

    log_div <- part.test.response*log_div
    if (progress) {utils::setTxtProgressBar(pb, 13)}

    part.test.response <- NULL

    deviances <- log_div - part.test.residuals
    if (progress) {utils::setTxtProgressBar(pb, 14)}

    part.test.residuals <- NULL; log_div <- NULL

    dev[fold,] <- 2*Matrix::colMeans(deviances)
    if (progress) {utils::setTxtProgressBar(pb, 15)}

    deviances <- NULL

  }

  s_min <- max(which(colMeans(dev) == min(colMeans(dev))))

  names <- list(sample_list = sample_list, gene_list = gene_list, mut_types_list = mut_types_list, col_names = col_names)
  return(list(fit = fit, dev = dev, s_min = s_min, names = names))

}

#' Fit Generative Model Without Sample-Specific Effects
#'
#' @description A function to fit a generative model to a mutation dataset that does not incorporate sample-specific effects.
#' Otherwise acts similarly to the function fit_gen_model().
#'
#' NOTE: fits produced by this model will not be compatible with predictive model fits downstream - it is purely for comparing
#' with full models.
#'
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param matrix (Matrix::sparseMatrix)
#' A mutation matrix, such as produced by the function get_table_from_maf().
#' @param sample_list (character)
#' The set of samples to be modelled.
#' @param gene_list (character)
#' The set of genes to be modelled.
#' @param mut_types_list (character)
#' The set of mutation types to be modelled.
#' @param col_names (character)
#' The column names of the 'matrix' parameter.
#' @param table (list)
#' Optional parameter combining matrix, sample_list, gene_list, mut_types_list, col_names, as is produced by the function get_tables().
#' @param nlambda (numeric)
#' The length of the vector of penalty weights, passed to the function glmnet::glmnet().
#' @param n_folds (numeric)
#' The number of cross-validation folds to employ.
#' @param maxit (numeric)
#' Technical parameter passed to the function glmnet::glmnet().
#' @param seed_id (numeric)
#' Input value for the function set.seed().
#' @param progress (logical)
#' Show progress bars and text.
#'
#' @return A list comprising three objects:
#' * An object 'fit', a fitted glmnet model.
#' * A table 'dev', giving average deviances for each regularisation penalty factor and cross-validation fold.
#' * An integer 's_min', the index of the regularsisation penalty minimising cross-validation deviance.
#' * A list 'names', containing the sample, gene, and mutation type information of the training data.
#' @export
#'
#' @examples
#' example_gen_model_unisamp <- fit_gen_model_unisamp(example_maf_data$gene_lengths,
#'                                                    table = example_tables$train)
#' print(names(example_gen_model))

fit_gen_model_unisamp <- function(gene_lengths, matrix = NULL, sample_list = NULL, gene_list = NULL, mut_types_list = NULL, col_names = NULL,
                                  table = NULL, nlambda = 100, n_folds = 10, maxit = 1e9, seed_id = 1234, progress = FALSE) {

  set.seed(seed_id)
  if (progress) {trace.it = 1}
  else {trace.it = 0}

  if (is.null(table) & any(is.null(matrix), is.null(sample_list), is.null(gene_list), is.null(mut_types_list), is.null(col_names))) {
    stop("If not providing a full tables object, must provide the inputs sample_list, gene_list, mut_types_list and col_names")
  }

  if(!is.null(table)) {
    matrix <- table$matrix
    sample_list <- table$sample_list
    gene_list <- table$gene_list
    mut_types_list <- table$mut_types_list
    col_names <- table$col_names
  }

  n_samples <- length(sample_list)
  n_genes <- length(gene_list)
  n_mut_types <- length(mut_types_list)

  if (any(dim(matrix) != c(n_samples, n_genes * n_mut_types))) {
    stop(paste0("Matrix has dimension ", dim(matrix), ", should have dimension ", c(n_samples, n_genes * n_mut_types)))
  }

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

  # Making output vector
  mutation_vector <- matrix
  dim(mutation_vector) <- c(n_samples * n_genes * n_mut_types, 1)

  # Constructing design matrix
  i <- 1:(n_samples*n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types), each = n_samples)
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_samples*n_mut_types)

  design_matrix <- Matrix::sparseMatrix(c(i, i), c(j, j_))

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  # We use weights (rather than offsets) to fit the model, as it makes glmnet more efficient.
  weights <- rep(t_s, each = n_samples, times = n_genes)
  weights <- weights*rep(gene_lengths[gene_list,]$max_cds, each = n_samples*n_mut_types)
  weights <- weights/mean(weights)

  weighted_observations <- ifelse(weights == 0, 0, as.vector(mutation_vector)/weights)

  # First (full dataset) run
  if (progress) {print("First glmnet run (full dataset)")}
  fit <- glmnet::glmnet(x = design_matrix, y = weighted_observations, nlambda = nlambda, weights = weights, family = "poisson", trace.it = trace.it, maxit = maxit)

  # Cross-validation
  if (progress) {print("Cross-validation:")}
  partitions <- sample(rep(1:n_folds, n_samples * n_genes * n_mut_types / n_folds), n_samples * n_genes * n_mut_types)
  dev <- matrix(0, n_folds, length(fit$lambda))

  for (fold in 1:n_folds) {
    part.design <- design_matrix[partitions != fold, ]
    part.response <- weighted_observations[partitions != fold]
    part.weights <- weights[partitions != fold]

    if (progress) {writeLines(paste("\nFitting glmnet on fold", fold))}
    part.fit <- glmnet::glmnet(part.design, y = part.response, family = "poisson", lambda = fit$lambda, weights = part.weights, trace.it = trace.it, maxit = maxit)
    part.weights <- NULL; part.design <- NULL; part.response <- NULL;

    if (progress) {print("Computing statistics")}
    ### This code should be made much simpler - it was originally written to squeeze the maximum amount out of a laptop without crashing.
    if (progress) {pb <- utils::txtProgressBar(max = 15, width = 100, style = 3)}

    part.test.design <- design_matrix[partitions == fold,]
    if (progress) {utils::setTxtProgressBar(pb, 1)}

    part.test.product <- part.test.design %*% part.fit$beta
    if (progress) {utils::setTxtProgressBar(pb, 2)}

    part.test.design <- NULL;

    part.test.product <- part.test.product + part.fit$a0
    if (progress) {utils::setTxtProgressBar(pb, 3)}

    part.fit <- NULL;

    part.test.product <- exp(part.test.product)
    if (progress) {utils::setTxtProgressBar(pb, 4)}

    part.test.weights <- Matrix::Matrix(rep(weights[partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda))
    if (progress) {utils::setTxtProgressBar(pb, 5)}

    part.test.predict <- part.test.weights * part.test.product
    if (progress) {utils::setTxtProgressBar(pb, 6)}

    part.test.product <- NULL;  part.test.weights <- NULL

    part.test.response <- Matrix::Matrix(rep(mutation_vector[1:(n_samples*n_mut_types*n_genes),][partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    if (progress) {utils::setTxtProgressBar(pb, 7)}

    part.test.residuals <- part.test.response - part.test.predict
    if (progress) {utils::setTxtProgressBar(pb, 8)}

    part.test.logpredict <- ifelse(part.test.predict == 0, 0, log(part.test.predict))
    if (progress) {utils::setTxtProgressBar(pb, 9)}

    part.test.predict <- NULL

    log_response <- mutation_vector[partitions == fold, ]
    log_response <- ifelse(log_response == 0, 0, log(log_response))
    if (progress) {utils::setTxtProgressBar(pb, 10)}

    part.test.logresponse <- Matrix::Matrix(rep(log_response, length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    if (progress) {utils::setTxtProgressBar(pb, 11)}

    log_response <- NULL

    log_div <- part.test.logresponse - part.test.logpredict
    if (progress) {utils::setTxtProgressBar(pb, 12)}

    part.test.logresponse <- NULL; part.test.logpredict <- NULL;

    log_div <- part.test.response*log_div
    if (progress) {utils::setTxtProgressBar(pb, 13)}

    part.test.response <- NULL

    deviances <- log_div - part.test.residuals
    if (progress) {utils::setTxtProgressBar(pb, 14)}

    part.test.residuals <- NULL; log_div <- NULL

    dev[fold,] <- 2*Matrix::colMeans(deviances)
    if (progress) {utils::setTxtProgressBar(pb, 15)}

    deviances <- NULL

  }

  s_min <- max(which(colMeans(dev) == min(colMeans(dev))))

  names <- list(sample_list = sample_list, gene_list = gene_list, mut_types_list = mut_types_list, col_names = col_names)
  return(list(fit = fit, dev = dev, s_min = s_min, names = names))

}



#' Fit Generative Model Without Gene/Variant Type-Specific Interactions
#'
#' @description A function to fit a generative model to a mutation dataset that does not incorporate gene/variant-specific effects.
#' Otherwise acts similarly to the function fit_gen_model().
#'
#' NOTE: fits produced by this model will not be compatible with predictive model fits downstream - it is purely for comparing
#' with full models.
#'
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param matrix (Matrix::sparseMatrix)
#' A mutation matrix, such as produced by the function get_table_from_maf().
#' @param sample_list (character)
#' The set of samples to be modelled.
#' @param gene_list (character)
#' The set of genes to be modelled.
#' @param mut_types_list (character)
#' The set of mutation types to be modelled.
#' @param col_names (character)
#' The column names of the 'matrix' parameter.
#' @param table (list)
#' Optional parameter combining matrix, sample_list, gene_list, mut_types_list, col_names, as is produced by the function get_tables().
#' @param nlambda (numeric)
#' The length of the vector of penalty weights, passed to the function glmnet::glmnet().
#' @param n_folds (numeric)
#' The number of cross-validation folds to employ.
#' @param maxit (numeric)
#' Technical parameter passed to the function glmnet::glmnet().
#' @param seed_id (numeric)
#' Input value for the function set.seed().
#' @param progress (logical)
#' Show progress bars and text.
#'
#' @return A list comprising three objects:
#' * An object 'fit', a fitted glmnet model.
#' * A table 'dev', giving average deviances for each regularisation penalty factor and cross-validation fold.
#' * An integer 's_min', the index of the regularsisation penalty minimising cross-validation deviance.
#' * A list 'names', containing the sample, gene, and mutation type information of the training data.
#' @export
#'
#' @examples
#' example_gen_model_unisamp <- fit_gen_model_unisamp(example_maf_data$gene_lengths,
#'                                                    table = example_tables$train)
#' print(names(example_gen_model))

fit_gen_model_uninteract <- function(gene_lengths, matrix = NULL, sample_list = NULL, gene_list = NULL, mut_types_list = NULL, col_names = NULL,
                                     table = NULL, nlambda = 100, n_folds = 10, maxit = 1e9, seed_id = 1234, progress = FALSE) {
  set.seed(seed_id)
  if (progress) {trace.it = 1}
  else {trace.it = 0}

  if (is.null(table) & any(is.null(matrix), is.null(sample_list), is.null(gene_list), is.null(mut_types_list), is.null(col_names))) {
    stop("If not providing a full tables object, must provide the inputs sample_list, gene_list, mut_types_list and col_names")
  }

  if(!is.null(table)) {
    matrix <- table$matrix
    sample_list <- table$sample_list
    gene_list <- table$gene_list
    mut_types_list <- table$mut_types_list
    col_names <- table$col_names
  }

  n_samples <- length(sample_list)
  n_genes <- length(gene_list)
  n_mut_types <- length(mut_types_list)

  if (any(dim(matrix) != c(n_samples, n_genes * n_mut_types))) {
    stop(paste0("Matrix has dimension ", dim(matrix), ", should have dimension ", c(n_samples, n_genes * n_mut_types)))
  }

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

  # Making output vector
  mutation_vector <- matrix
  dim(mutation_vector) <- c(n_samples * n_genes * n_mut_types, 1)

  # Constructing design matrix
  i <- 1:(n_samples*n_genes*n_mut_types)
  j <- rep(1:n_genes, each = n_samples * n_mut_types)

  design_matrix <- Matrix::sparseMatrix(i, j)

  t_i_getter <- Matrix::sparseMatrix(i = rep(1:n_samples, n_genes*n_mut_types), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_samples, n_samples*n_genes*n_mut_types))
  t_i <- as.vector(t_i_getter %*% mutation_vector)
  t_i_getter <- NULL

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  # We use weights (rather than offsets) to fit the model, as it makes glmnet more efficient.
  weights <- rep(t_s, each = n_samples, times = n_genes)* rep(t_i, times = n_mut_types*n_genes)
  weights <- weights*rep(gene_lengths[gene_list,]$max_cds, each = n_samples*n_mut_types)
  weights <- weights/mean(weights)

  weighted_observations <- ifelse(weights == 0, 0, as.vector(mutation_vector)/weights)

  # First (full dataset) run
  if (progress) {print("First glmnet run (full dataset)")}
  fit <- glmnet::glmnet(x = design_matrix, y = weighted_observations, nlambda = nlambda, weights = weights, family = "poisson", trace.it = trace.it, maxit = maxit)

  # Cross-validation
  if (progress) {print("Cross-validation:")}
  partitions <- sample(rep(1:n_folds, n_samples * n_genes * n_mut_types / n_folds), n_samples * n_genes * n_mut_types)
  dev <- matrix(0, n_folds, length(fit$lambda))

  for (fold in 1:n_folds) {
    part.design <- design_matrix[partitions != fold, ]
    part.response <- weighted_observations[partitions != fold]
    part.weights <- weights[partitions != fold]

    if (progress) {writeLines(paste("\nFitting glmnet on fold", fold))}
    part.fit <- glmnet::glmnet(part.design, y = part.response, family = "poisson", lambda = fit$lambda, weights = part.weights, trace.it = trace.it, maxit = maxit)
    part.weights <- NULL; part.design <- NULL; part.response <- NULL;

    if (progress) {print("Computing statistics")}
    ### This code should be made much simpler - it was originally written to squeeze the maximum amount out of a laptop without crashing.
    if (progress) {pb <- utils::txtProgressBar(max = 15, width = 100, style = 3)}

    part.test.design <- design_matrix[partitions == fold,]
    if (progress) {utils::setTxtProgressBar(pb, 1)}

    part.test.product <- part.test.design %*% part.fit$beta
    if (progress) {utils::setTxtProgressBar(pb, 2)}

    part.test.design <- NULL;

    part.test.product <- part.test.product + part.fit$a0
    if (progress) {utils::setTxtProgressBar(pb, 3)}

    part.fit <- NULL;

    part.test.product <- exp(part.test.product)
    if (progress) {utils::setTxtProgressBar(pb, 4)}

    part.test.weights <- Matrix::Matrix(rep(weights[partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda))
    if (progress) {utils::setTxtProgressBar(pb, 5)}

    part.test.predict <- part.test.weights * part.test.product
    if (progress) {utils::setTxtProgressBar(pb, 6)}

    part.test.product <- NULL;  part.test.weights <- NULL

    part.test.response <- Matrix::Matrix(rep(mutation_vector[1:(n_samples*n_mut_types*n_genes),][partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    if (progress) {utils::setTxtProgressBar(pb, 7)}

    part.test.residuals <- part.test.response - part.test.predict
    if (progress) {utils::setTxtProgressBar(pb, 8)}

    part.test.logpredict <- ifelse(part.test.predict == 0, 0, log(part.test.predict))
    if (progress) {utils::setTxtProgressBar(pb, 9)}

    part.test.predict <- NULL

    log_response <- mutation_vector[partitions == fold, ]
    log_response <- ifelse(log_response == 0, 0, log(log_response))
    if (progress) {utils::setTxtProgressBar(pb, 10)}

    part.test.logresponse <- Matrix::Matrix(rep(log_response, length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    if (progress) {utils::setTxtProgressBar(pb, 11)}

    log_response <- NULL

    log_div <- part.test.logresponse - part.test.logpredict
    if (progress) {utils::setTxtProgressBar(pb, 12)}

    part.test.logresponse <- NULL; part.test.logpredict <- NULL;

    log_div <- part.test.response*log_div
    if (progress) {utils::setTxtProgressBar(pb, 13)}

    part.test.response <- NULL

    deviances <- log_div - part.test.residuals
    if (progress) {utils::setTxtProgressBar(pb, 14)}

    part.test.residuals <- NULL; log_div <- NULL

    dev[fold,] <- 2*Matrix::colMeans(deviances)
    if (progress) {utils::setTxtProgressBar(pb, 15)}

    deviances <- NULL

  }

  s_min <- max(which(colMeans(dev) == min(colMeans(dev))))

  names <- list(sample_list = sample_list, gene_list = gene_list, mut_types_list = mut_types_list, col_names = col_names)
  return(list(fit = fit, dev = dev, s_min = s_min, names = names))

}


#' Investigate Generative Model Comparisons
#'
#' @description Given a generative model of the type we propose, and an alternate version (saturated "S", sample-independent "US", gene-independent "UG" or gene/variant interaction independent "UI"), either produces the estimated observations on the training dataset or calculates residual deviance between models.
#'
#' @param training_data (list) Likely the 'train' component of a call to get_mutation_tables().
#' @param gen_model (list) A generative model - result of a call to fit_gen_model*().
#' @param alt_gen_model (list) An alternative generative model.
#' @param alt_model_type (character) One of "S" (saturated), "US" (sample-independent), "UG", (gene-independent), "UI" (gene/variant-interaction independent).
#' @param gene_lengths (dataframe) A gene lengths data frame.
#' @param calculate_deviance (logical) If TRUE, returns residual deviance statistics. If FALSE, returns training data predictions.
#'
#' @return If calculate_deviance = FALSE:
#'
#' A list with two entries, est_mut_vec and alt_est_mut_vec, each of length n_samples x n_genes x n_mut_types, giving expected mutation value for each combination of sample, gene and variant type in the training dataset under the two models being compared.
#'
#' If calculate_deviance = TRUE:
#'
#' A list with two entries, deviance and df, corresponding to the residual deviance and residual degrees of freedom between the two models on the training set.
#' @export
#'
#' @examples sat_dev <- get_gen_estimates(training_data = example_tables$train,
#'                                        gen_model = example_gen_model,
#'                                        alt_model_type = "S",
#'                                        gene_lengths = example_maf_data$gene_lengths,
#'                                        calculate_deviance = TRUE)

get_gen_estimates <- function(training_data, gen_model, alt_gen_model = NULL, alt_model_type = "S", gene_lengths = NULL, calculate_deviance = FALSE) {

  n_samples <- length(training_data$sample_list)
  n_genes <- length(training_data$gene_list)
  n_mut_types <- length(training_data$mut_types_list)

  mutation_vector <- training_data$matrix
  dim(mutation_vector) <- c(n_samples * n_genes * n_mut_types, 1)
  mutation_vector <- as.vector(mutation_vector)

  t_i_getter <- Matrix::sparseMatrix(i = rep(1:n_samples, n_genes*n_mut_types), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_samples, n_samples*n_genes*n_mut_types))
  t_i <- as.vector(t_i_getter %*% mutation_vector)
  t_i_getter <- NULL

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  # We use weights (rather than offsets) to fit the model, as it makes glmnet more efficient.
  weights <- rep(t_s, each = n_samples, times = n_genes)* rep(t_i, times = n_mut_types*n_genes)
  weights <- weights*rep(gene_lengths[training_data$gene_list,]$max_cds, each = n_samples*n_mut_types)
  weights <- weights/mean(weights)

  i <- 1:(n_samples*n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types), each = n_samples)
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_samples*n_mut_types)

  design_matrix <- Matrix::sparseMatrix(c(i, i), c(j, j_))

  est_mut_vec <- exp(gen_model$fit$a0[gen_model$s_min] + design_matrix %*% gen_model$fit$beta[, gen_model$s_min])
  est_mut_vec <- as.vector(est_mut_vec * weights)

  rm(i); rm(j); rm(j_);

 if (alt_model_type == "US") {
   weights_us <- rep(t_s, each = n_samples, times = n_genes)
   weights_us <- weights_us*rep(gene_lengths[training_data$gene_list,]$max_cds, each = n_samples*n_mut_types)
   weights_us <- weights_us/mean(weights_us)

   alt_est_mut_vec <- exp(gen_model$fit$a0[alt_gen_model$s_min] + design_matrix %*% alt_gen_model$fit$beta[, alt_gen_model$s_min])
   alt_est_mut_vec <- as.vector(alt_est_mut_vec * weights_us)
   }

  else if (alt_model_type == "UG") {
    alt_est_mut_vec <- rep(gene_lengths[training_data$gene_list,]$max_cds, each = n_samples*n_mut_types)
    alt_est_mut_vec <- alt_est_mut_vec * rep(t_s, each = n_samples, times = n_genes)
    alt_est_mut_vec <- alt_est_mut_vec / (n_mut_types * n_genes * mean(alt_est_mut_vec))
    alt_est_mut_vec <- alt_est_mut_vec * rep(rep(t_i, times = n_mut_types*n_genes))
  }

  else if (alt_model_type == "UI") {
  i <- 1:(n_samples*n_genes*n_mut_types)
  j <- rep(1:n_genes, each = n_samples * n_mut_types)
  design_matrix_ui <- Matrix::sparseMatrix(i,j)

  alt_est_mut_vec <- exp(gen_model$fit$a0[alt_gen_model$s_min] + design_matrix_ui %*% alt_gen_model$fit$beta[, alt_gen_model$s_min])
  alt_est_mut_vec <- as.vector(alt_est_mut_vec * weights)
  }

  else if (alt_model_type %in% c("S")){

  }
  else {
    stop(paste(alt_model_type, "not a valid alternative model type."))
  }

  if (calculate_deviance) {
    deviance <- 2*sum(mutation_vector*ifelse(mutation_vector == 0, 0, log(mutation_vector/est_mut_vec)) - (mutation_vector - est_mut_vec))
    if (alt_model_type == "S") {
      df <- n_samples * n_genes * n_mut_types - (n_samples +  n_mut_types + sum(gen_model$fit$beta[, gen_model$s_min] != 0))
      return(list(deviance = deviance, df = df))
    }
    else {
      alt_deviance <- 2*sum(mutation_vector*ifelse(mutation_vector == 0, 0, log(mutation_vector/alt_est_mut_vec)) - (mutation_vector - alt_est_mut_vec))

      if (alt_model_type == "US") {
        df <- n_samples
      }

      else if (alt_model_type == "UG") {
        df <- sum(gen_model$fit$beta[, gen_model$s_min] != 0)
      }

      else if (alt_model_type == "UI") {
        df <- sum(gen_model$fit$beta[, gen_model$s_min] != 0) - sum(alt_gen_model$fit$beta[, alt_gen_model$s_min] != 0)
      }
    }
    return(list(deviance = alt_deviance - deviance, df = df))
  }
  else {
    return(list(est_mut_vec = est_mut_vec, alt_est_mut_vec = alt_est_mut_vec))
  }
}


#' Visualise Generative Model Fit
#'
#' @description A function to visualise how well a general model has fitted to a mutation dataset across cross-validation folds. Designed to produce a similar output to glmnet's function plot.cv.glmnet.
#'
#' @param gen_model (list)
#' A generative model fitted by fit_gen_model()
#' @param x_sparsity
#' Show model sparsity on x axis rather than lambda.
#' @param y_sparsity
#' Show model sparsity on y axis rather than deviance.
#' @param mut_type
#' Produce separate plots for each mutation type.
#'
#' @return
#' Summary plot of the generative model fit across folds.
#'
#' @export
#'
#' @examples
#' p <- vis_model_fit(example_gen_model)
#'

vis_model_fit <- function(gen_model, x_sparsity = FALSE, y_sparsity = FALSE, mut_type = NULL) {
  fit_data <- data.frame(log_lambda = log(gen_model$fit$lambda),
                         nzero = Matrix::colSums(gen_model$fit$beta != 0),
                         Deviance = colMeans(gen_model$dev),
                         sd = matrixStats::colSds(gen_model$dev))

  fit_data$ymin <- fit_data$Deviance - fit_data$sd
  fit_data$ymax <- fit_data$Deviance + fit_data$sd
  fit_data$colour <- (fit_data$Deviance == min(fit_data$Deviance))

  if (!x_sparsity & !y_sparsity) {
    p <- ggplot2::ggplot(fit_data, ggplot2::aes_string(x = "log_lambda", y = "Deviance", ymin = "ymin",
                                                ymax = "ymax", colour = "colour")) +
      ggplot2::geom_pointrange() + ggplot2::labs(x = latex2exp::TeX("$\\log(\\kappa_1)$"), y = "Average Deviance") +
      ggplot2::scale_colour_manual(values = c("black", "red")) +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none") + ggplot2::theme(axis.title.y = ggplot2::element_text(vjust = 2))
  }
  else {

  }

  return(p)
}
