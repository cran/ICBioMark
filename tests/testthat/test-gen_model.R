context("Generative Model Fitting")

test_that("A simple call to fit_gen_model() has the right format and value of output", {
  output <- fit_gen_model(gene_lengths = example_maf_data$gene_lengths, table = example_tables$train)

  expect_is(output, "list")
  expect_equal(names(output), c("fit", "dev", "s_min", "names"))
})


