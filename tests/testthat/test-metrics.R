context("Functions for analysing the output of predictive models")


test_that("get_r_squared gives the right format of output", {
  output <- get_r_squared(predictions = example_predictions,
                          biomarker_values = example_tmb_tables$val,
                          model = "Refitted T")

  expect_type(output, "list")
  expect_equal(colnames(output), c("panel_length", "model", "biomarker", "stat", "metric"))
})

test_that("get_auprc gives the right format of output", {
  output <- get_auprc(predictions = example_predictions,
                          biomarker_values = example_tmb_tables$val,
                          model = "Refitted T")

  expect_type(output, "list")
  expect_equal(colnames(output), c("panel_length", "model", "biomarker", "stat", "metric"))
})

test_that("get_stats gives the right format of output", {
  output <- get_stats(predictions = example_predictions,
                      biomarker_values = example_tmb_tables$val,
                      model = "Refitted T")

  expect_type(output, "list")
  expect_equal(colnames(output), c("panel_length", "model", "biomarker", "stat", "metric"))
})
