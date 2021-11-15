#' Construct Optimisation Parameters.
#'
#' @description An internal function. From the learned generative model and training data, produces a vector of weights p to be used in
#' the subsequent group lasso optimisation, alongside a biomarker-dependent normalisation quantity p_norm.
#'
#' @param gen_model (list)
#' A generative mutation model, fitted by fit_gen_model().
#' @param training_matrix (sparse matrix)
#' A sparse matrix of mutations in the training dataset, produced by get_mutation_tables().
#' @param marker_mut_types (character)
#' A character vector listing which mutation types (of the set specified in the generative model
#' attribute 'names') constitute the biomarker in question.
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#'
#' @return
#' A list with three entries:
#'  * A vector p, with an entry corresponding to each combination of gene and mutation type specified
#'   in the generative model fitted. Each component is a non-negative value corresponding to a weighting
#'   p to be supplied to a group lasso optimisation.
#'  * A numeric p_norm, giving the factor between p_{gs} and phi_{0gs} (see paper for details).
#'  * A vector biomarker_columns, detailing which of the elements of p correspond to gene/mutation type
#'   combinations contributing to the biomarker in question.
#'
#' @export
#'
#' @examples
#' p <- get_p(example_gen_model, example_tables$train$matrix,
#'            marker_mut_types = c("I"), gene_lengths = example_maf_data$gene_lengths)
#' print(p$p[1:5])
#' print(p$p_norm)
#' print(p$bc[1:5])


get_p <- function(gen_model, training_matrix, marker_mut_types, gene_lengths) {

  n_samples <- nrow(training_matrix)
  n_genes <- length(gen_model$names$gene_list)
  n_mut_types <- length(gen_model$names$mut_types_list)

  mutation_vector <- training_matrix
  dim(mutation_vector) <- c(n_samples*n_genes*n_mut_types, 1)
  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

  biomarker_ids <- purrr::map(marker_mut_types, ~which(gen_model$names$mut_types_list == .))
  biomarker_columns <- sort(unlist(purrr::map(biomarker_ids, ~seq(., n_genes*n_mut_types, n_mut_types))))

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples),
                                     j = 1:(n_samples*n_genes*n_mut_types),
                                     dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  weights <- rep(t_s, times = n_genes)
  weights <- weights*rep(gene_lengths[gen_model$names$gene_list,]$max_cds, each = n_mut_types)
  weights <- weights/mean(weights)

  i <- 1:(n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types))
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_mut_types)

  param_getter <- Matrix::sparseMatrix(c(i, i), c(j, j_))
  param_est <- as.vector(param_getter %*% gen_model$fit$beta[, gen_model$s_min])
  names(param_est) <- gen_model$names$col_names
  param_getter <- NULL

  p_norm <- sum((weights * exp(gen_model$fit$a0[gen_model$s_min] + param_est))[biomarker_columns])
  p <- weights * exp(gen_model$fit$a0[gen_model$s_min] + param_est) / p_norm
  names(p) <- gen_model$names$col_names

  return(list(p = p, p_norm = p_norm, bc = biomarker_columns))
}

#' Construct Bias Penalisation
#'
#' @description An internal function, producing the correct bias penalisation for use in predictive model fitting.
#'
#' @param gen_model (list)
#' A generative mutation model, fitted by fit_gen_model().
#' @param p_norm (numeric)
#' Scaling factor between coefficients of p and parameters of generative model (see paper for details).
#' @param training_matrix (sparse matrix)
#' A sparse matrix of mutations in the training dataset, produced by get_mutation_tables().
#' @param marker_training_values (dataframe)
#' A dataframe containing training values for the biomarker in question.
#' @param method (function)
#' How to select a representative biomarker value from the training dataset. Defaults to max().
#'
#' @return
#' A numerical value, to be used as a penalty weighting in the subsequent group lasso optimisation.
#' @export
#'
#' @examples
#' K <- get_K(example_gen_model, 1, example_tables$train$matrix)
#' print(K)

get_K <- function(gen_model, p_norm, training_matrix, marker_training_values = NULL, method = max) {

  if (is.null(marker_training_values)) {
    marker_training_values <- data.frame(Tumor_Sample_Barcode = gen_model$names$sample_list,
                                         value = Matrix::rowSums(training_matrix),
                                         stringsAsFactors = FALSE)
  }

  K <- exp(gen_model$fit$a0[[gen_model$s_min]])*method(marker_training_values$value)*p_norm
  return(K)

}

#' Extract Panel Details from Group Lasso Fit
#'
#' @description An internal function for analysing a group Lasso fit as part of the predictive model learning procedure, which returns the sets of genes identified by different iterations of the group Lasso algorithm.
#'
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param fit (list)
#' A fit from the group lasso algorithm, produced by the function gglasso (package: gglasso).
#' @param gene_list (character)
#' A character vector of genes listing the genes (in order) included in the model pred_fit.
#' @param mut_types_list (character)
#' A character vector listing the mutation type groupings (in order) included in the model pred_fit.
#'
#' @return
#' A list of two elements:
#'  * panel_genes: A matrix where each row corresponds to a gene, each column to an iteration of the group
#'   lasso with a different penalty factor, and the elements booleans specifying whether that gene was selected
#'   to be included in that iteration.
#'  * panel_lengths:
#' @export
#'
#' @examples
#' panels <- get_panels_from_fit(example_maf_data$gene_lengths, example_first_pred_tmb$fit,
#' example_gen_model$names$gene_list, mut_types_list = example_gen_model$names$mut_types_list)
#'
#' print(panels$fit)

get_panels_from_fit <- function(gene_lengths, fit, gene_list, mut_types_list) {

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol
  n_genes <- length(gene_list)
  n_mut_types <- length(mut_types_list)
  col_names <- paste0(rep(gene_list, each = n_mut_types), "_", rep(mut_types_list, times = n_genes))
  pred_fit_abs <- as.data.frame(abs(fit$beta))
  panel_inclusion <- dplyr::group_by(pred_fit_abs, gene = factor(rep(gene_list,
                                            each = length(mut_types_list)),
                                        levels = gene_list))
  panel_genes <- dplyr::select(dplyr::summarise(panel_inclusion, dplyr::across(1:length(fit$lambda), sum)), setdiff(colnames(panel_inclusion), "gene")) != 0
  rownames(panel_genes) <- gene_list
  panel_gene_lengths <- gene_lengths[gene_list,]$max_cds

  panel_lengths <- as.vector(t(panel_gene_lengths) %*% panel_genes)

  return(list(panel_genes = panel_genes, panel_lengths = panel_lengths))

}

#' First-Fit Predicitve Model with Group Lasso
#'
#' @description This function implements the first-fit procedure described in Bradley and Cannings, 2021. It requires at least a generative model and a dataframe containing gene lengths as input.
#'
#' @param gen_model (list)
#' A generative mutation model, fitted by fit_gen_model().
#' @param lambda (numeric)
#' A vector of penalisation weights for input to the group lasso optimiser gglasso.
#' @param biomarker (character)
#' The biomarker in question. If "TMB" or "TIB", then automatically defines the subsequent
#' variable marker_mut_types.
#' @param marker_mut_types (character)
#' The set of mutation type groupings constituting the biomarker being estimated. Should be
#' a vector comprising of elements of the mut_types_list vector in the 'names' attribute of
#' gen_model.
#' @param training_matrix (sparse matrix)
#' A sparse matrix of mutations in the training dataset, produced by get_mutation_tables().
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param marker_training_values (dataframe)
#' A dataframe containing two columns: 'Tumor_Sample_Barcode', containing the sample IDs for the training
#' dataset, and a second column containing training values for the biomarker in question.
#' @param K_method (function)
#' How to select a representative biomarker value from the training dataset. Defaults to max().
#' @param free_genes (character)
#' Which genes should escape penalisation (for example when augmenting a pre-existing panel).
#'
#' @return
#' A list of six elements:
#'  * fit: Output of call to gglasso.
#'  * panel_genes: A matrix where each row corresponds to a gene, each column to an iteration of the group
#'   lasso with a different penalty factor, and the elements booleans specifying whether that gene was selected
#'   to be included in that iteration.
#'  * panel_lengths: A vector giving total panel length for each gglasso iteration.
#'  * p: The vector of weights used in the optimisation procedure.
#'  * K: The bias penalty factor used in the optimisation procedure.
#'  * names: Gene and mutation type information as used when fitting the generative model.
#'
#' @export
#'
#' @examples
#' example_first_fit <- pred_first_fit(example_gen_model, lambda = exp(seq(-9, -14, length.out = 100)),
#'                                     training_matrix = example_tables$train$matrix,
#'                                     gene_lengths = example_maf_data$gene_lengths)
#'
pred_first_fit <- function(gen_model, lambda = exp(seq(-16,-24, length.out = 100)), biomarker = "TMB",
                           marker_mut_types = c("NS", "I"), training_matrix, gene_lengths, marker_training_values = NULL,
                           K_method = max, free_genes = c()) {

  if (biomarker == "TIB") {
    marker_mut_types <- c("I")
  }

  wrong_mutation_types <- setdiff(marker_mut_types, gen_model$names$mut_types_list)
  if (length(wrong_mutation_types) > 0) {
    stop(paste0("Mutation types ", paste(wrong_mutation_types, collapse = ", "), " not in generative model."))
  }

  n_samples <- nrow(training_matrix)
  n_genes <- length(gen_model$names$gene_list)
  n_mut_types <- length(gen_model$names$mut_types_list)

  message("Getting p")
  p <- get_p(gen_model = gen_model, training_matrix = training_matrix, marker_mut_types = marker_mut_types,
             gene_lengths = gene_lengths)
  message("Getting K")
  K <- get_K(gen_model = gen_model, p_norm = p$p_norm, training_matrix = training_matrix,
             marker_training_values = marker_training_values, method = K_method)

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol
  reduced_gene_lengths <- gene_lengths
  reduced_gene_lengths[intersect(free_genes, gene_lengths$Hugo_Symbol),'max_cds'] <- 0
  pf <- reduced_gene_lengths[gen_model$names$gene_list,]$max_cds
  message("Making matrix")
  X <- matrix(0, n_mut_types * n_genes + 1, n_mut_types * n_genes)
  colnames(X)

  X[1,] <- sqrt(K)*p$p
  diag(X[2:(n_mut_types * n_genes + 1),]) <- sqrt(p$p)
  Y <- c(sqrt(K), sqrt(p$p))
  Y[setdiff(2:(n_mut_types * n_genes + 1), 1 + p$bc)] <- 0

  message("Fitting")
  fit <- gglasso::gglasso(x = X, y = Y, loss = "ls", lambda = lambda,
                          group = rep(1:n_genes, each = n_mut_types), intercept = FALSE, pf = pf)
  rownames(fit$beta) <- gen_model$names$col_names

  panels <- get_panels_from_fit(gene_lengths = gene_lengths, fit = fit,
                                gene_list = gen_model$names$gene_list, mut_types_list = gen_model$names$mut_types_list)

  return(list(fit = fit, panel_genes = panels$panel_genes, panel_lengths = panels$panel_lengths, p = p$p, K = K, names = gen_model$names))
}

#' Refitted Predictive Model for a Given Panel
#'
#' @description A function taking the output of a call to pred_first_fit(), as well as gene length information, and a specified panel (list of genes), and producing a refitted predictive model on that given panel.
#'
#' @param pred_first (list)
#' A first-fit predictive model as produced by pred_first_fit().
#' @param gene_lengths (dataframe)
#' A dataframe of gene lengths (see example_maf_data$gene_lengths for format).
#' @param model (character)
#' A choice of "T", "OLM" or "Count" specifying how predictions should be made.
#' @param genes (character)
#' A vector of gene names detailing the panel being used.
#' @param biomarker (character)
#' If "TMB" or "TIB", automatically defines marker_mut_types, otherwise this will
#' need to be specified separately.
#' @param marker_mut_types (character)
#' A vector specifying which mutation types groups determine the biomarker in question.
#' @param training_data (list)
#' Training data, as produced by get_mutation_tables() (select train, val or test).
#' @param training_values (dataframe)
#' Training true values, as produced by get_biomarker_tables() (select train, val or test).
#' @param mutation_vector (numeric)
#' Optional vector specifying the values of the training matrix (training_data$matrix) in
#' vector rather than matrix form.
#' @param t_s (numeric)
#' Optional vector specifying the frequencies of different mutation types.
#'
#' @return
#' A list with three elements:
#'  * fit, a list including a sparse matrix 'beta' giving prediction weights.
#'  * panel_genes, a sparse (logical) matrix giving the genes included in prediction.
#'  * panel_lengths, a singleton vector giving the length of the panel used.
#' @export
#'
#' @examples
#' example_refit_panel <- pred_refit_panel(pred_first = example_first_pred_tmb,
#'   gene_lengths = example_maf_data$gene_lengths, genes = paste0("GENE_", 1:10))

pred_refit_panel <- function(pred_first = NULL, gene_lengths = NULL, model = "T", genes, biomarker = "TMB",
                             marker_mut_types = c("NS", "I"), training_data = NULL,
                             training_values = NULL, mutation_vector = NULL, t_s = NULL) {
  n_genes <- length(pred_first$names$gene_list)

  if (!is.null(training_data)) {
    n_mut_types <- length(training_data$mut_types)
  }

  else {
    n_mut_types <- length(pred_first$names$mut_types_list)
  }

  wrong_genes_model <- setdiff(genes, pred_first$names$gene_list)
  if (length(wrong_genes_model) > 0) {
    warning(paste("Eliminating the following genes not in the generative model: ", paste0(wrong_genes_model, collapse = ", ")))
    genes <- intersect(genes, pred_first$names$gene_list)
    }

  wrong_genes_lengths <- setdiff(genes, gene_lengths$Hugo_Symbol)
  if(length(wrong_genes_lengths > 0)) {
    warning(paste("Eliminating the following genes without gene lengths: ", paste0(wrong_genes_lengths, collapse = ", ")))
    genes <- intersect(genes, gene_lengths)
  }

  if (biomarker == "TIB") {
    marker_mut_types <- c("I")
  }
  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol
  if (length(genes) > 0) {
    panel_lengths <- c(sum(gene_lengths[genes,]$max_cds))
  }
  else {
    panel_lengths <- c(0)
  }

  if (!is.null(training_data)) {
    cols_panel <- paste0(rep(genes, each = n_mut_types), "_", training_data$mut_types_list)
  }
  else{
    cols_panel <- paste0(rep(genes, each = n_mut_types), "_", pred_first$names$mut_types_list)
  }

  panel_genes <- Matrix::Matrix(pred_first$names$col_names %in% cols_panel,
                                nrow = n_genes * n_mut_types, ncol = 1,
                                sparse = TRUE)

  if (model == "T") {
    if (is.null(pred_first)) {
      stop("Need first-fit model (pred_first) for fitting T estimator.")
    }

    if (length(genes) > 0) {
      p_panel <- pred_first$p[cols_panel]
      p_reduced <- unlist(purrr::map(1:n_mut_types, ~sum(p_panel[seq(., length(p_panel), n_mut_types)])))
      names(p_reduced) <- pred_first$names$mut_types_list
      X_panel <- diag(x = p_reduced, nrow = length(p_reduced)) + pred_first$K * p_reduced %*% t(p_reduced)
      Y_panel <- (pred_first$K + (pred_first$names$mut_types_list %in% marker_mut_types))*p_reduced
      w_panel <- solve(X_panel) %*% Y_panel
      beta <- Matrix::Matrix(0, nrow = n_genes * n_mut_types, ncol = 1, sparse = TRUE)
      rownames(beta) <- pred_first$names$col_names
      beta[cols_panel, ] <- rep(w_panel, times = length(genes))
      fit <- list(beta = beta)
    }
    else {
      beta <- Matrix::Matrix(0, nrow = n_genes * n_mut_types, ncol = 1, sparse = TRUE)
      fit <- list(beta = beta)
    }
  }
  else if (model == "OLM") {
    n_samples <- nrow(training_data$matrix)
    if (is.null(training_data) | is.null(training_values)) {
      stop("Need training matrix and values for OLM fitting.")
    }
    colnames(training_data$matrix) <- training_data$col_names
    rownames(training_data$matrix) <- training_values$Tumor_Sample_Barcode
    training_values <- training_values[pred_first$names$sample_list,]
    train <- as.data.frame(cbind(matrix(training_values[[biomarker]], n_samples, 1),
                   as.matrix(training_data$matrix[,cols_panel])))
    colnames(train) <- c(biomarker, cols_panel)
    formula <- stats::as.formula(paste(biomarker, "~ . - 1"))
    fit <- stats::lm(formula = formula, data = train)

    fit$beta <- Matrix::Matrix(0, n_genes * n_mut_types, sparse = TRUE)
    rownames(fit$beta) <- training_data$col_names
    fit$beta[cols_panel,] <- fit$coefficients
    fit$beta[is.na(fit$beta)] <- 0

  }
  else if (model == "Count") {
    if (is.null(training_data) | is.null(training_values)) {
      stop("Need training matrix and values for Count fitting.")
    }

    n_samples <- nrow(training_data$matrix)
    n_genes <- length(pred_first$names$gene_list)
    n_mut_types <- length(training_data$mut_types_list)

    rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

    if (is.null(mutation_vector)) {
      mutation_vector <- training_data$matrix
      dim(mutation_vector) <- c(n_samples*n_genes*n_mut_types, 1)
    }

    if (is.null(t_s)) {
      t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples),
                                         j = 1:(n_samples*n_genes*n_mut_types),
                                         dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
      t_s <- as.vector(t_s_getter %*% mutation_vector)
      t_s_getter <- NULL
    }

    lengths_factor <- sum(gene_lengths[pred_first$names$gene_list, 'max_cds']) / sum(gene_lengths[genes, 'max_cds'])
    mut_types_factor <- sum(t_s[pred_first$names$mut_types_list %in% marker_mut_types])/ sum(t_s)

    beta <- Matrix::Matrix(0, nrow = n_genes * n_mut_types, ncol = 1, sparse = TRUE)
    rownames(beta) <- training_data$col_names
    beta[cols_panel,] <- lengths_factor * mut_types_factor

    fit <- list(beta = beta)
  }

  else {
    stop("Model should be one of 'T', 'OLM', or 'Count'")
  }
  return(list(fit = fit, panel_genes = panel_genes, panel_lengths = panel_lengths))
}

#' Get Refitted Predictive Models for a First-Fit Range of Panels
#'
#' @description  A function producing a refitted predictive model for each panel produced by usage of the function pred_first_fit(), by repeatedly applying the function pred_refit_panel().
#'
#' @param pred_first (list)
#' A first-fit predictive model as produced by pred_first_fit().
#' @param gene_lengths (dataframe)
#' A dataframe of gene lengths (see example_maf_data$gene_lengths for format).
#' @param model (character)
#' A choice of "T", "OLM" or "Count" specifying how predictions should be made.
#' @param biomarker (character)
#' If "TMB" or "TIB", automatically defines marker_mut_types, otherwise this will
#' need to be specified separately.
#' @param marker_mut_types (character)
#' A vector specifying which mutation types groups determine the biomarker in question.
#' @param training_data (sparse matrix)
#' Training matrix, as produced by get_mutation_tables() (select train, val or test).
#' @param training_values (dataframe)
#' Training true values, as produced by get_biomarker_tables() (select train, val or test).
#' @param mutation_vector (numeric)
#' Optional vector specifying the values of the training matrix (training_data$matrix) in
#' vector rather than matrix form.
#' @param t_s (numeric)
#' Optional vector specifying the frequencies of different mutation types.
#' @param max_panel_length (numeric)
#' Upper bound for panels to fit refitted models to. Most useful for "OLM" and "Count"
#' model types.
#'
#' @return
#'  A list with three elements:
#'  * fit, a list including a sparse matrix 'beta' giving prediction weights for each
#'  first-fit panel (one panel per column).
#'  * panel_genes, a sparse (logical) matrix giving the genes included in prediction
#'  for each first-fit panel.
#'  * panel_lengths, a vector giving the length of each first-fit panel.
#' @export
#'
#' @examples
#' example_refit_range <- pred_refit_range(pred_first = example_first_pred_tmb,
#'   gene_lengths = example_maf_data$gene_lengths)
pred_refit_range <- function(pred_first = NULL, gene_lengths = NULL, model = "T", biomarker = "TMB",
                             marker_mut_types = c("NS", "I"), training_data = NULL, training_values = NULL,
                             mutation_vector = NULL, t_s = NULL, max_panel_length = NULL) {

    if (model == "Count") {

      if (is.null(training_data) | is.null(training_values)) {
        stop("Need training matrix and values for Count fitting.")
      }
      n_mut_types <- length(training_data$mut_types_list)
      n_samples <- nrow(training_data$matrix)
      n_genes <- length(training_data$gene_list)

      if (is.null(mutation_vector)) {
        mutation_vector <- training_data$matrix
        dim(mutation_vector) <- c(n_samples*n_genes*n_mut_types, 1)
      }

      if (is.null(t_s)) {
        t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples),
                                           j = 1:(n_samples*n_genes*n_mut_types),
                                           dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
        t_s <- as.vector(t_s_getter %*% mutation_vector)
        t_s_getter <- NULL
      }
    }
    which_genes <- which(pred_first$panel_genes, arr.ind = TRUE)
    genes <- purrr::map(1:ncol(pred_first$panel_genes), ~ rownames(which_genes[which_genes[, 'col'] == ., ]))

    if (!is.null(max_panel_length)) {
      s <- max(which(pred_first$panel_lengths <= max_panel_length))
      genes <- genes[1:s]
      pred_first$panel_genes <- pred_first$panel_genes[,1:s, drop = FALSE]
      pred_first$panel_lengths <- pred_first$panel_lengths[1:s]
    }

    betas <- purrr::map(genes, ~ pred_refit_panel(genes = ., pred_first = pred_first, gene_lengths = gene_lengths, model = model,
                                                  biomarker = biomarker, marker_mut_types = marker_mut_types,
                                                  training_data = training_data, training_values = training_values,
                                                  mutation_vector = mutation_vector, t_s = t_s)$fit$beta)

    beta <- do.call(cbind, betas)
  return(list(fit = list(beta = beta), panel_genes = pred_first$panel_genes, panel_lengths = pred_first$panel_lengths))

}

#' Produce Predictions on an Unseen Dataset
#'
#' @description A function taking a predictive model(s) and new observations, and applying the predictive model to them to return predicted biomarker values.
#'
#' @param pred_model (list)
#' A predictive model as fitted by pred_first_fit(), pred_refit_panel() or
#' pred_refit_range().
#' @param new_data (list)
#' A new dataset, containing a matrix of observations and a list of sample IDs.
#' Likely comes from the 'train', 'val' or 'test' argument of a call to
#' get_mutation_tables().
#' @param s (numeric)
#' If producing predictions for a single panel, s chooses which panel
#' (column in a pred_fit object) to produce predictions for.
#' @param max_panel_length (numeric)
#' If producing predictions for a single panel, maximum panel length to
#' specify that panel.
#'
#' @return
#' A list with two elements:
#' * predictions, a matrix containing a row for each sample and a column for each
#' panel.
#' * panel_lengths, a vector containing the length of each panel.
#' @export
#'
#' @examples
#' example_predictions <- get_predictions(example_refit_range, new_data =
#' example_tables$val)

get_predictions <- function(pred_model, new_data,
                            s = NULL, max_panel_length = NULL) {
  predictions <- as.matrix(new_data$matrix %*% pred_model$fit$beta)
  rownames(predictions) <- new_data$sample_list

  if (!is.null(max_panel_length)) {
    s = max(which(pred_model$panel_lengths <= max_panel_length))
  }
  if (!is.null(s)) {
    predictions <- predictions[, s, drop = FALSE]
  }
  return(list(predictions = predictions, panel_lengths = pred_model$panel_lengths))
}

#' Produce Error Bounds for Predictions
#'
#' @description A function to produce a confidence region for a linear predictor. In upcoming versions will (hopefully) be greatly simplified.
#'
#' @param predictions (list)
#' A predictions object, as produced by get_predictions().
#' @param pred_model (list)
#' A predictive model, as produced by pred_first_fit(), pred_refit_panel() or
#' pred_refit_range().
#' @param gen_model (list)
#' A generative model, as produce by fit_gen_model
#' @param training_matrix (sparse matrix)
#' A training matrix, as produced by get_tables()$matrix or get_table_from_maf()$matrix.
#' @param gene_lengths (data frame)
#' A data frame with columns 'Hugo_Symbol' and 'max_cds'. See example_maf_data$gene_lengths, or
#' ensembl_gene_lengths for examples.
#' @param biomarker_values (data frame)
#' A data frame containing the true values of the biomarker in question.
#' @param alpha (numeric)
#' Confidence level for error bounds.
#' @param range_factor (numeric)
#' Value specifying how far beyond the range of max(biomarker) to plot confidence region.
#' @param s (numeric)
#' If input predictions are for a range of panels, s chooses which panel
#' (column in a pred_fit object) to produce predictions for.
#' @param max_panel_length (numeric)
#' Select panel by maximum length.
#' @param biomarker (character)
#' Which biomarker is being predicted.
#' @param marker_mut_types (character)
#' If biomarker is not one of "TMB" or "TIB", then this is required to specify which mutation type
#' groups constitute the biomarker.
#' @param model (character)
#' The model (must be based on a linear estimator) for which prediction intervals are being generated.
#'
#' @return
#' A list with two entries:
#'  * prediction_intervals:
#'  * confidence_region:
#' @export
#'
#' @examples
#' example_intervals <- pred_intervals(predictions = get_predictions(example_refit_range,
#'                new_data = example_tables$val),
#'                pred_model = example_refit_range, biomarker_values = example_tmb_tables$val,
#'                gen_model = example_gen_model, training_matrix = example_tables$train$matrix,
#'                max_panel_length = 15000, gene_lengths = example_maf_data$gene_lengths)
#'
#' example_confidence_plot <- ggplot2::ggplot() +
#'   ggplot2::geom_point(data = example_intervals$prediction_intervals,
#'              ggplot2::aes(x = true_value, y = estimated_value)) +
#'         ggplot2::geom_ribbon(data = example_intervals$confidence_region,
#'           ggplot2::aes(x = x, ymin = y_lower, ymax = y_upper),
#'                     fill = "red", alpha = 0.2) +
#'         ggplot2::geom_line(data = example_intervals$confidence_region,
#'           ggplot2::aes(x = x, y = y), linetype = 2) +
#'         ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
#'
#' plot(example_confidence_plot)

pred_intervals <- function(predictions, pred_model, gen_model, training_matrix, gene_lengths,
                           biomarker_values, alpha = 0.1, range_factor = 1.1, s = NULL, max_panel_length = NULL,
                           biomarker = "TMB", marker_mut_types = c("NS", "I"),
                           model = "Refitted T") {
  if (!is.null(max_panel_length)) {
    s = max(which(predictions$panel_lengths <= max_panel_length))
  }
  if (!is.null(s)) {
    predictions$predictions <- predictions$predictions[, s, drop = FALSE]
  }


  else if (ncol(predictions$predictions > 1)) {
    stop("Either provide a predictions matrix with one column, or select which column
         should be used with s or max_panel_length.")
  }

  if (ncol(pred_model$fit$beta) > 1) {
    warning(paste0("Using column ",  s, " of predictive model fit"))
  }

  if (biomarker == "TIB") {
    marker_mut_types = c("I")
  }

  pred_range <- seq(min(biomarker_values[[biomarker]]) , range_factor * max(biomarker_values[[biomarker]]), length.out = 100)

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

  n_samples <- nrow(training_matrix)
  n_genes <- length(gen_model$names$gene_list)
  n_mut_types <- length(gen_model$names$mut_types_list)

  mutation_vector <- training_matrix
  dim(mutation_vector) <- c(n_samples*n_genes*n_mut_types, 1)

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples),
                                     j = 1:(n_samples*n_genes*n_mut_types),
                                     dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  i <- 1:(n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types))
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_mut_types)

  param_getter <- Matrix::sparseMatrix(c(i, i), c(j, j_))
  param_est <- as.vector(param_getter %*% gen_model$fit$beta[, gen_model$s_min])
  names(param_est) <- gen_model$names$col_names
  param_getter <- NULL

  exp_g_s <- exp(param_est)*rep(gene_lengths[gen_model$names$gene_list, 'max_cds'], each = n_mut_types)*rep(t_s, n_genes)
  exp_mu <- Matrix::rowSums(training_matrix) / (exp(gen_model$fit$a0[gen_model$s_min]) * sum(exp_g_s))

  biomarker_ids <- purrr::map(marker_mut_types, ~which(gen_model$names$mut_types_list == .))
  biomarker_columns <- sort(unlist(purrr::map(biomarker_ids, ~seq(., n_genes*n_mut_types, n_mut_types))))

  bias <- (sum(exp_g_s[biomarker_columns]) - sum(exp_g_s*as.vector(pred_model$fit$beta[,s])))^2
  var <- sum(exp_g_s[biomarker_columns]*(1 - as.vector(pred_model$fit$beta[,s])[biomarker_columns])^2) +
    sum(exp_g_s[-biomarker_columns]*as.vector(pred_model$fit$beta[,s])[-biomarker_columns]^2)
  norm <- sum(exp_g_s*as.vector(pred_model$fit$beta[,s]))

  prediction_intervals <- data.frame(Tumor_Sample_Barcode = biomarker_values$Tumor_Sample_Barcode,
                                 true_value = biomarker_values[[biomarker]],
                                 estimated_value = predictions$predictions[biomarker_values$Tumor_Sample_Barcode, ],
                                 model = model)
  prediction_intervals$upper <-  (2*prediction_intervals$true_value + var/(alpha*norm) +
                                    sqrt((2*prediction_intervals$true_value + var/(alpha*norm))^2 - 4*(1- bias/(alpha*norm^2))*prediction_intervals$true_value^2)) /
                                    (2*(1 - bias / (alpha*norm^2)))
  prediction_intervals$lower <- (2*prediction_intervals$true_value + var/(alpha*norm) -
                                   sqrt((2*prediction_intervals$true_value + var/(alpha*norm))^2 - 4*(1- bias/(alpha*norm^2))*prediction_intervals$true_value^2)) /
                                   (2*(1 - bias / (alpha*norm^2)))


  confidence_region <- data.frame(x = pred_range, y = pred_range,
                              y_lower = (2*pred_range + var/(alpha*norm) -
                                           sqrt((2*pred_range + var/(alpha*norm))^2 - 4*(1- bias/(alpha*norm^2))*pred_range^2)) /
                                           (2*(1- bias / (alpha*norm^2))),
                              y_upper = (2*pred_range + var/(alpha*norm) +
                                           sqrt((2*pred_range + var/(alpha*norm))^2 - 4*(1- bias/(alpha*norm^2))*pred_range^2)) /
                                           (2*(1- bias / (alpha*norm^2))),
                              model = model)

  return(list(prediction_intervals = prediction_intervals, confidence_region = confidence_region))
}
