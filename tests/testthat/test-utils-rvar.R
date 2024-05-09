test_that("rvar utils work", {
  S <- 2000
  x <- posterior::rvar(rnorm(S, mean = 1, sd = 1))
  x0 <- rvar_set(x, 2.3)
  x1 <- rvar_std_normal(x)
  expect_equal(mean(x0)[1], 2.3)
  expect_true(all(abs(mean(x1)) < 0.1))

  n <- 4
  mu <- rep(1:n, each = S)
  y <- posterior::rvar(array(rnorm(S * n, mean = mu, sd = 2), dim = c(S, n)))
  y0 <- rvar_set(y, -5.4)
  y1 <- rvar_std_normal(y)
  expect_true(all(mean(y0) == -5.4))
  expect_true(all(abs(mean(y1)) < 0.1))
})
