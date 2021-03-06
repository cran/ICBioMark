% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{get_auprc}
\alias{get_auprc}
\title{AUPRC Metrics for Predictions}
\usage{
get_auprc(predictions, biomarker_values, model = "", threshold = 300)
}
\arguments{
\item{predictions}{(list)
A list with two elements, 'predictions' and 'panel_lengths',
as produced by the function get_predictions().}

\item{biomarker_values}{(dataframe)
A dataframe with two columns, 'Tumor_Sample_Barcode' and a
column with the name of the biomarker in question containing values.}

\item{model}{(character)
The name of the model type producing these predictions.}

\item{threshold}{(numeric)
The threshold for biomarker high/low categorisation.}
}
\value{
A dataframe with 5 columns:
\itemize{
\item panel_length: the length of each panel.
\item model: the model that produced the predictions.
\item biomarker: the name of the biomarker in question.
\item stat: the AUPRC values for each panel.
\item metric: a constant character "AUPRC".
}
}
\description{
A function to return AUPRC metrics for predictions vs actual values. Works well when piped to straight from get_predictions().
}
\examples{
example_auprc <- get_auprc(predictions = get_predictions(example_refit_panel,
new_data = example_tables$val), biomarker_values = example_tmb_tables$val,
model = "Refitted T", threshold = 10)
}
