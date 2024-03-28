test_that(".bf_weighted_sum_coloc works", {
  cur_out <- .bf_weighted_sum_coloc(w = rep(0.25,4),
                                    bf = seq(1,4),
                                    fp_coloc = 0.5,
                                    i = 1)

  expect_equal(length(cur_out),1)
})
