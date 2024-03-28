test_that("multi_intact works", {
  data("multi_simdat")

  cur_out <- multi_intact(df = multi_simdat)

  expect_equal(nrow(cur_out[[1]]),nrow(multi_simdat))
})
