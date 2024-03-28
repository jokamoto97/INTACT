test_that(".multi_prior_estimation works", {
  data("multi_simdat")
  multi_simdat$fp_coloc <- linear(pmax(multi_simdat$GLCP_1,multi_simdat$GLCP_2))

  pi0 <- 1 - .pi1_fun(z_vec = qnorm(pchisq(multi_simdat$chisq,
                                           df = 2,lower.tail = FALSE)/2))
  pi_start <- c(pi0,rep(1-pi0,3)/3)

  cur_out <- .multi_prior_estimation(df = multi_simdat,
                                     pi_init = pi_start,
                                     chisq_vec = multi_simdat$chisq,
                                     z_1 = multi_simdat$z_1,
                                     z_2 = multi_simdat$z_2,
                                     fp_coloc = multi_simdat$fp_coloc)

  expect_equal(nrow(cur_out[[1]]),nrow(multi_simdat))
})
