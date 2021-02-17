context("Inputting MAFs, outputting mutation matrices.")

test_that("A simple call to get_table_from_maf() has the right structure and value", {
  output <- get_table_from_maf(maf = example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))

  expect_is(output, "list")
  expect_equal(length(names(output)), 5)
  expect_equal(names(output), c("matrix", "sample_list", "gene_list", "mut_types_list", "col_names"))
  expect_known_output(output, file = "datatest_table.rds")
})

test_that("A simple call to get_mutation_tables() has the right structure and value", {
  output <- get_mutation_tables(maf = example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100), gene_list = paste0("GENE_", 1:20))

  expect_is(output, "list")
  expect_equal(length(names(output)), 3)
  expect_equal(names(output), c("train", "val", "test"))

  train_data <- output$train
  expect_is(train_data, "list")
  expect_equal(length(names(train_data)), 5)
  expect_equal(names(train_data), c("matrix", "sample_list", "gene_list", "mut_types_list", "col_names"))

  expect_is(train_data$matrix, "dgCMatrix")

  expect_equal(output, example_tables)
})

test_that("Simple calls to get_biomarker_tables() has the right structure and value", {
  tmb <- get_biomarker_tables(maf = example_maf_data$maf, biomarker = "TMB", sample_list = paste0("SAMPLE_", 1:100))
  tib <- get_biomarker_tables(maf = example_maf_data$maf, biomarker = "TIB", sample_list = paste0("SAMPLE_", 1:100))

  expect_is(tmb, "list")
  expect_is(tib, "list")

  expect_equal(names(tmb), c("train", "val", "test"))
  expect_equal(names(tib), c("train", "val", "test"))

  expect_equal(colnames(tmb$train), c("Tumor_Sample_Barcode", "TMB"))
  expect_equal(colnames(tib$train), c("Tumor_Sample_Barcode", "TIB"))

  expect_equal(tmb, example_tmb_tables)
  expect_equal(tib, example_tib_tables)
})
