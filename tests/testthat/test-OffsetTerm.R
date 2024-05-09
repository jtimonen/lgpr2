test_that("different types of OffsetTerm work correctly", {
  m <- LonModel$new(y ~ offset(id2 | arm) + gp(time))
  dat <- testdata
  dat$id2 <- dat$id
  fit <- m$fit(dat, chains = 1)
  expect_s3_class(fit$plot(), "ggplot")
})
