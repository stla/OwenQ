context("spowen2")

test_that("Tolerance factors", {
  n <- 10
  ke <- 3.197
  t <- ke*sqrt(n)
  expect_equal(spowen2(n-1, t, delta), 0.95, tolerance=1e-3)
  ke <- 2.84
  t <- ke*sqrt(n)
  expect_equal(spowen2(n-1, t, delta), 0.9, tolerance=1e-3)
  ke <- 4.066
  t <- ke*sqrt(n)
  expect_equal(spowen2(n-1, t, delta), 0.99, tolerance=1e-3)
})
