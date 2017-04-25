context("OwenQ2")

test_that("OwenQ2 for R=0 equals ptOwen", {
  owenQ2 <- sapply(1:8, function(nu) OwenQ2(nu=nu, t=2, delta=2, R=0))
  owenS <- sapply(1:8, function(nu) ptOwen(q=2, nu=nu, delta=2))
  expect_equal(owenQ2, owenS, tolerance=1e-15)
})

test_that("OwenQ1+OwenQ2 equals ptOwen", {
  owenQ2 <- sapply(1:8,
                   function(nu){
                     OwenQ1(nu=nu, t=2, delta=2, R=1) + OwenQ2(nu=nu, t=2, delta=2, R=1)
                     })
  owenS <- sapply(1:8, function(nu) ptOwen(q=2, nu=nu, delta=2))
  expect_equal(owenQ2, owenS, tolerance=1e-15)
})

test_that("OwenQ2 for t=+Inf does not depend on delta", {
  expect_true(OwenQ2(5, Inf, 2, 2) == OwenQ2(5, Inf, 3, 2))
  expect_true(OwenQ2(6, Inf, 2, 2) == OwenQ2(6, Inf, 3, 2))
  # does not depend on t for delta=-Inf and the same result
  expect_true(OwenQ2(5, Inf, 2, 2) == OwenQ2(5, 1, -100, 2))
  expect_true(OwenQ2(6, Inf, 2, 2) == OwenQ2(6, 1, -100, 2))
})

test_that("OwenQ2 for t=-Inf equals 0", {
  expect_true(OwenQ2(5, -Inf, 2, 1) == 0)
  expect_true(OwenQ2(6, -Inf, 2, 1) == 0)
})

test_that("OwenQ2 for |t|=Inf and R=0", {
  expect_true(is.nan(OwenQ2(1, Inf, 2, 0))) # should give 1
  expect_true(is.nan(OwenQ2(1, -Inf, 2, 0))) # should give 0
})


test_that("OwenQ2 - bivariate Student", {
  t1 <- 2; t2 <- 1; delta1 <- 3; delta2 <- 2
  nu <- 6
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  owen <- - OwenQ2(nu, t2, delta2, R) + OwenQ2(nu, t1, delta1, R)
  wolfram <- 0.03257737810540227
  expect_equal(owen, wolfram, tolerance=1e-10)
  #
  nu <- 5
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  owen <- - OwenQ2(nu, t2, delta2, R) + OwenQ2(nu, t1, delta1, R)
  wolfram <- 0.0353568969628651
  expect_equal(owen, wolfram, tolerance=1e-9)
})

