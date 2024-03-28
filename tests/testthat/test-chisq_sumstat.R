test_that("chisq_sumstat works", {
  data(z_sumstats)
  data(protwt_sumstats)
  data(exprwt_sumstats)
  data(ld_sumstats)

  cur_out <- chisq_sumstat(z_vec = z_sumstats,
                           w = cbind(protwt_sumstats,exprwt_sumstats),
                           R = ld_sumstats)
  expect_equal(length(cur_out),1)
})
