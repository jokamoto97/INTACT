test_that(".logistic_em works", {
  data(simdat)
  cur_out <- .logistic_em(d_vec = sample(c(0,1),1197,replace=TRUE),
                          pprobs = intact(GLCP_vec=simdat$GLCP, prior_fun=linear,
                          z_vec = simdat$TWAS_z, t = 0.05),
                          alpha = c(0,0))
  expect_equal(length(cur_out),2)
})
