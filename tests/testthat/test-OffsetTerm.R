test_that("different types of OffsetTerm work correctly", {
  m <- TSModel$new(y ~ offset(id2 | arm) + gp(t))
  dat <- testdata
  dat$id2 <- dat$id
  fit <- m$fit(dat, chains = 1)
  expect_s3_class(fit$plot(), "ggplot")
})
