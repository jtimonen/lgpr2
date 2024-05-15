test_that("projection predictive search works", {
  dat <- simulate(40, 0.3, c(0, 1, 0, 0), 0.2)
  df <- data.frame(cbind(dat$x, dat$y))
  colnames(df)[ncol(df)] <- "y"

  mod <- LonModel$new(y ~ gp(x1) + gp(x2) + gp(x3) + gp(x4))
  fit <- mod$fit(data = df, chains = 1)

  fs <- pp_forward_search(fit)

  plt1 <- plot_pp_pexp(fs)
  plt2 <- plot_pp_elpd(fs, fit$loo_estimate())
  expect_s3_class(plt1, "ggplot")
  expect_s3_class(plt2, "ggplot")
})
