test_priors <- function() {
  checkEquals(step(0.2,t = 0.5), 0)
  checkEquals(hybrid(1,1))
  checkEquals(linear(0.5,0.5))
  checkEquals(expit(0.03,0))
}
