test_that(".enrich_bootstrap_se works", {
  data(simdat)
  cur_out <- .enrich_bootstrap_se(d_vec = sample(c(0,1),1197,replace=TRUE),
                                  pprobs = intact(GLCP_vec=simdat$GLCP,
                                                  prior_fun=linear,
                                                  z_vec = simdat$TWAS_z,
                                                  t = 0.05))
  expect_equal(sum(cur_out[3,] == 1),100)
})
