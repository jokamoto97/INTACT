test_that(".multi_em_posteriors", {
  cur_out <- .multi_em_posteriors(w = rep(0.25,4),
                                    bf = seq(1,4),
                                    fp_coloc = 0.5)

  expect_equal(length(cur_out),4)
})
