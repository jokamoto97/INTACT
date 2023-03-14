test_that("fdr_rst works", {
  data(simdat)
  cur_out <- fdr_rst(simdat$GLCP)

  expect_equal(nrow(cur_out), nrow(simdat))
})
