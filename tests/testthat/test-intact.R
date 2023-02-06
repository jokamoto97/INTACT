test_that("intact works", {
  data(simdat)
  cur_out <- intact(GLCP_vec=simdat$GLCP, z_vec = simdat$TWAS_z)

  expect_equal(sum(cur_out == 1),18)
})
