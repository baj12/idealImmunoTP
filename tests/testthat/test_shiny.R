library("idealImmunoTP")

context("Check that shiny app is generated")

dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
# rlt <- rlogTransformation(dds)
# cm <- counts(dds)
# cd <- colData(dds)

test_that("Shiny app is generated", {
  expect_is(idealImmunoTP(), "shiny.appobj")
  expect_is(idealImmunoTP(dds_obj = dds), "shiny.appobj")
})
