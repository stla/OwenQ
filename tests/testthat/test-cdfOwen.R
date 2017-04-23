context("cdfOwen")

test_that("pOwen4", {
  t1 <- 2; t2 <- 1; delta1 <- 3; delta2 <- 2
  nu <- 6
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  diff <- OwenQ1(nu, t2, delta2, R) - OwenQ1(nu, t1, delta1, R)
  owen4 <- pOwen4(nu, t1, t2, delta1, delta2)
  expect_equal(diff, owen4, tolerance=1e-17)
  wolfram <- 0.01785518085912236
  expect_equal(owen4, wolfram, tolerance=1e-11)
  #
  nu <- 5
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  diff <- OwenQ1(nu, t2, delta2, R) - OwenQ1(nu, t1, delta1, R)
  owen4 <- pOwen4(nu, t1, t2, delta1, delta2)
  expect_true(diff == owen4)
  wolfram <- 0.01868982415809893
  expect_equal(owen4, wolfram, tolerance=1e-9)
})
