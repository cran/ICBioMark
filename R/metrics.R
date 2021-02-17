#' R Squared Metrics for Predictions
#'
#' @description A function to return R^2 metrics for predictions vs actual values. Works well when piped to straight from get_predictions().
#'
#' @param predictions (list)
#' A list with two elements, 'predictions' and 'panel_lengths',
#' as produced by the function get_predictions().
#' @param biomarker_values (dataframe)
#' A dataframe with two columns, 'Tumor_Sample_Barcode' and a
#' column with the name of the biomarker in question containing values.
#' @param model (character)
#' The name of the model type producing these predictions.
#' @param threshold (numeric)
#' Unusued in this function: present for calls to get_stats().
#'
#' @return
#' A dataframe with 5 columns:
#'  * panel_length: the length of each panel.
#'  * model: the model that produced the predictions.
#'  * biomarker: the name of the biomarker in question.
#'  * stat: the R squared values for each panel.
#'  * metric: a constant character "R" for R squared.
#' @export
#'
#' @examples
#' example_r <- get_r_squared(predictions = get_predictions(example_refit_panel, new_data =
#'   example_tables$val), biomarker_values = example_tmb_tables$val, model = "Refitted T")

get_r_squared <- function(predictions, biomarker_values, model = "", threshold = 10) {

  predictions$predictions <- predictions$predictions[biomarker_values$Tumor_Sample_Barcode, , drop = FALSE]

  if (length(predictions$panel_lengths) != ncol(predictions$predictions)) {
    stop("panel_lengths doesn't match predictions")
  }

  if (nrow(predictions$predictions) != nrow(biomarker_values)) {
    stop("predictions doesn't match biomarker_values")
  }

  n_panels <- length(predictions$panel_lengths)
  biomarker <- colnames(biomarker_values)[2]

  r.squared <- purrr::map(1:n_panels, ~ 1 - sum((predictions$predictions[, .] - biomarker_values[[biomarker]])^2) /
                sum((biomarker_values[[biomarker]] - mean(biomarker_values[[biomarker]]))^2))
  output <- data.frame(panel_length = predictions$panel_lengths,
                       model = model, biomarker = biomarker, stat = unlist(r.squared), metric = "R")
  return(output)

  }

#' AUPRC Metrics for Predictions
#'
#' @description A function to return AUPRC metrics for predictions vs actual values. Works well when piped to straight from get_predictions().
#'
#' @param predictions (list)
#' A list with two elements, 'predictions' and 'panel_lengths',
#' as produced by the function get_predictions().
#' @param biomarker_values (dataframe)
#' A dataframe with two columns, 'Tumor_Sample_Barcode' and a
#' column with the name of the biomarker in question containing values.
#' @param model (character)
#' The name of the model type producing these predictions.
#' @param threshold (numeric)
#' The threshold for biomarker high/low categorisation.
#'
#' @return
#' A dataframe with 5 columns:
#'  * panel_length: the length of each panel.
#'  * model: the model that produced the predictions.
#'  * biomarker: the name of the biomarker in question.
#'  * stat: the AUPRC values for each panel.
#'  * metric: a constant character "AUPRC".
#' @export
#'
#' @examples
#' example_auprc <- get_auprc(predictions = get_predictions(example_refit_panel,
#' new_data = example_tables$val), biomarker_values = example_tmb_tables$val,
#' model = "Refitted T", threshold = 10)

get_auprc <- function(predictions, biomarker_values, model = "", threshold = 300) {

  predictions$predictions <- predictions$predictions[biomarker_values$Tumor_Sample_Barcode, , drop = FALSE]

  if (length(predictions$panel_lengths) != ncol(predictions$predictions)) {
    stop("panel_lengths doesn't match predictions")
  }

  if (nrow(predictions$predictions) != nrow(biomarker_values)) {
    stop("predictions doesn't match biomarker_values")
  }

  n_panels <- length(predictions$panel_lengths)
  biomarker <- colnames(biomarker_values)[2]
  classes <- biomarker_values[[biomarker]] >= threshold

  auprc <- purrr::map(1:n_panels, ~ PRROC::pr.curve(scores.class0 = predictions$predictions[classes, .],
                                             scores.class1 = predictions$predictions[!classes, .])$auc.integral)

  output <- data.frame(panel_length = predictions$panel_lengths,
                       model = model, biomarker = biomarker, stat = unlist(auprc), metric = "AUPRC")

  return(output)
}

#' Metrics for Predictive Performance
#'
#' @description A function to return a variety metrics for predictions vs actual values. Works well when piped to straight from get_predictions().
#'
#' @param predictions (list)
#' A list with two elements, 'predictions' and 'panel_lengths',
#' as produced by the function get_predictions().
#' @param biomarker_values (dataframe)
#' A dataframe with two columns, 'Tumor_Sample_Barcode' and a
#' column with the name of the biomarker in question containing values.
#' @param model (character)
#' The name of the model type producing these predictions.
#' @param threshold (numeric)
#' The threshold for biomarker high/low categorisation.
#' @param metrics (character)
#' A vector of the names of metrics to calculate.
#'
#' @return
#' dataframe with 5 columns:
#'  * panel_length: the length of each panel.
#'  * model: the model that produced the predictions.
#'  * biomarker: the name of the biomarker in question.
#'  * stat: the metric values for each panel.
#'  * metric: the name of the metric.
#' @export
#'
#' @examples
#' example_stat <- get_stats(predictions = get_predictions(example_refit_panel,
#' new_data = example_tables$val), biomarker_values = example_tmb_tables$val,
#' model = "Refitted T", threshold = 10)

get_stats <- function(predictions, biomarker_values, model = "", threshold = 300, metrics = c("R", "AUPRC")) {
  stats_list <- list(rep(NA, length(metrics)))
  index <- 1
  if ("R" %in% metrics) {
    stats_list[[index]] <- get_r_squared(predictions = predictions, biomarker_values = biomarker_values,
                                         model = model)
    index <- index + 1
  }

  if ("AUPRC" %in% metrics) {
    stats_list[[index]] <- get_auprc(predictions = predictions, biomarker_values = biomarker_values,
                                   model = model, threshold = threshold)
  }

  output <- do.call(rbind, stats_list)
  return(output)
}
