test_that("IdentityTransform works correctly", {
  x_test <- stats::rnorm(5)
  r <- IdentityTransform$new()
  a <- r$forward(x_test)
  expect_true(all(a == x_test))
  b <- r$backward(x_test)
  c <- r$backward(a)
  expect_true(all(b == x_test))
  expect_true(all(c == x_test))
})

test_that("LinearTransform works correctly", {
  x_test <- stats::rnorm(5)
  m <- rnorm(1)
  sd <- 3.4
  r <- LinearTransform$new(multiplier = sd, offset = m)
  expect_output(r$print(), "LinearTransform")
  a <- r$forward(x_test)
  b <- r$backward(a)
  expect_lt(max(b - x_test), 1e-12)
})


test_that("UnitScaleTransform works correctly", {
  x_test <- stats::rnorm(5)
  r <- UnitScaleTransform$new()
  r <- r$set_using_data(x_test)
  a <- r$forward(x_test)
  expect_output(r$print(), "UnitScaleTransform")
  expect_equal(max(a), 1.0)
  expect_equal(min(a), -1.0)
})

test_that("MaxScaleTransform works correctly", {
  x_test <- 100 * (2 + stats::rnorm(100))
  r <- MaxScaleTransform$new()
  r <- r$set_using_data(x_test)
  a <- r$forward(x_test)
  expect_output(r$print(), "MaxScaleTransform")
  expect_equal(max(a), 1.0)
})
