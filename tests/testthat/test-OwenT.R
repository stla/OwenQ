context("OwenT")

test_that("Owen T - reflection properties", {
  expect_true(OwenT(2,1) == OwenT(-2,1))
  expect_true(OwenT(2,1) == -OwenT(2,-1))
  expect_true(OwenT(2,0.1) == OwenT(-2,0.1))
  expect_true(OwenT(2,0.1) == -OwenT(2,-0.1))
  expect_true(OwenT(2,11) == OwenT(-2,11))
  expect_true(OwenT(2,11) == -OwenT(2,-11))
})

test_that("OwenT(0,a)", {
  a <- 1
  expect_true(OwenT(0,a) == atan(a)/(2*pi))
  a <- 0.9
  expect_true(OwenT(0,a) == atan(a)/(2*pi))
  a <- 2
  expect_true(OwenT(0,a) == atan(a)/(2*pi))
})

test_that("OwenT(h,1)", {
  h <- 2
  expect_equal(OwenT(h,1), pnorm(h)*(1-pnorm(h))/2, tolerance=1e-16)
  h <- 100 # >cut.point
  expect_true(OwenT(h,1) == pnorm(h)*(1-pnorm(h))/2) # this is 0
})

test_that("OwenT(h,Inf)", {
  h <- 1
  expect_equal(OwenT(h,Inf), (1-pnorm(abs(h)))/2, tolerance=1e-16)
  h <- 0.5
  expect_equal(OwenT(h,Inf), (1-pnorm(abs(h)))/2, tolerance=1e-16)
  h <- 10
  expect_equal(OwenT(h,Inf), (1-pnorm(abs(h)))/2, tolerance=1e-16)
  h <- 100
  expect_equal(OwenT(h,Inf), (1-pnorm(abs(h)))/2, tolerance=1e-16)
})

test_that("OwenT(Inf,a) = 0", {
  a <- 30
  expect_true(OwenT(Inf,a) == 0) # NaN
})

test_that("Relation OwenT Cauchy", {
  h <- 2; a <- 2
  expect_equal(OwenT(h, a), 1/2*(pt(a, 1, h*sqrt(1+a^2)) - pnorm(-h)),
               tolerance=1e-12)
})


test_that("Comparison Mathematica", {
  expect_equal(OwenT(1,2), 0.078468186993084, tolerance=1e-16)
  expect_equal(OwenT(1,0.5), 0.0430646911207853656324, tolerance=1e-17)
  expect_equal(OwenT(2,0.9), 0.0109285988291624569525, tolerance=1e-17)
})

test_that("Relation T(h,a) T(ah,1/a)", {
  h <- 0.5; a <- 2
  expect_equal(OwenT(h,a)+OwenT(a*h,1/a),
                (pnorm(h)+pnorm(a*h))/2-pnorm(h)*pnorm(a*h),
               tolerance=1e-16)
  a <- -2
  expect_equal(OwenT(h,a)+OwenT(a*h,1/a),
               (pnorm(h)+pnorm(a*h))/2-pnorm(h)*pnorm(a*h)-0.5,
               tolerance=1e-16)
  h <- 100; a <- 2 # this gives 0
  expect_true(OwenT(h,a)+OwenT(a*h,1/a) ==
                (pnorm(h)+pnorm(a*h))/2-pnorm(h)*pnorm(a*h))
})
