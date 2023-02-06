test_that("intactGSE works", {
  data(simdat)
  data(gene_set_list)
  cur_out <- intactGSE(gene_data = simdat,gene_sets = gene_set_list)

  expect_equal(sum(cur_out$CONVERGED),2)
})
