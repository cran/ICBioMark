context("Functions for fitting predictive models")

test_that("pred_first_fit gives reference value", {
  output <- pred_first_fit(example_gen_model, lambda = exp(seq(-9, -14, length.out = 100)),
                           training_matrix = example_tables$train$matrix,
                           gene_lengths = example_maf_data$gene_lengths)
  expect_equal(output, example_first_pred_tmb)
})

  test_that("pred_refit_panel gives the right format of output", {
    output <- pred_refit_panel(pred_first = example_first_pred_tmb,
                               gene_lengths = example_maf_data$gene_lengths,
                               genes = paste0("GENE_", 1:10))
    expect_type(output, "list")
    expect_length(output, 3)
    expect_named(output, c("fit", "panel_genes", "panel_lengths"))
    expect_type(output$fit, "list")
    expect_length(dim(output$panel_genes), 2)
    expect_true(typeof(output$panel_lengths) %in% c("integer", "double"))

    expect_equal(output, example_refit_panel)

})

  test_that("pred_refit_range gives the right format of output", {
    output <- pred_refit_range(pred_first = example_first_pred_tmb,
                               gene_lengths = example_maf_data$gene_lengths)

    expect_type(output, "list")
    expect_length(output, 3)
    expect_named(output, c("fit", "panel_genes", "panel_lengths"))
    expect_type(output$fit, "list")
    expect_length(dim(output$panel_genes), 2)
    expect_true(typeof(output$panel_lengths) %in% c("integer", "double"))

    expect_equal(output, example_refit_range)
})

  test_that("get_predictions gives the right format of output", {
    output <- get_predictions(example_refit_range, new_data =
                              example_tables$val)

    expect_type(output, "list")
    expect_length(output, 2)
    expect_named(output, c("predictions", "panel_lengths"))
    expect_length(dim(output$predictions), 2)
    expect_true(typeof(output$panel_lengths) %in% c("integer", "double"))
})
