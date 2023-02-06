test_that(".pi1_fun works", {
  data(simdat)
  cur_out <- .pi1_fun(simdat$TWAS_z)

  expect_equal(sum(is.numeric(cur_out)),1)
})
