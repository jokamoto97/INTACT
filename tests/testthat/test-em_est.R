test_that(".em_est works", {
  data(simdat)
  cur_out <- .em_est(d_vec = sample(c(0,1),1197,replace=TRUE),
                     pprobs = intact(GLCP_vec=simdat$GLCP,
                     prior_fun=linear, z_vec = simdat$TWAS_z, t = 0.05))

  expect_equal(sum(cur_out == 1),1)
})
