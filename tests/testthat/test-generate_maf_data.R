context("Simulating MAF data")

test_that("Check that the output of generate_maf_data() has the right structure", {
  output <- generate_maf_data()
  expect_is(output, "list")
  expect_equal(length(output), 2)
  expect_equal(names(output), c("maf", "gene_lengths"))
})
