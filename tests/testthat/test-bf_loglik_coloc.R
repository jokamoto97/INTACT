test_that(".bf_loglik_coloc", {
  cur_out <- .bf_loglik_coloc(w = rep(0.25,4),
                                    bf = seq(1,4),
                                    fp_coloc = 0.5)

  expect_equal(length(cur_out),1)
})
