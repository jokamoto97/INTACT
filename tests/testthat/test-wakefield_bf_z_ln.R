test_that(".wakefield_bf_z_ln works", {
  data(multi_simdat)
  cur_out <- .wakefield_bf_z_ln(z_vec = multi_simdat$z_1)

  expect_equal(length(cur_out),nrow(multi_simdat))
})
