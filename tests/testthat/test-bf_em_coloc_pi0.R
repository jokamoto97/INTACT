test_that(".bf_em_coloc_pi0 works", {

  cur_out <- .bf_em_coloc_pi0(w = rep(0.25,4),
                                    bf = rnorm(400),
                                    fp_coloc = runif(100))

  expect_equal(length(cur_out),4)
})
