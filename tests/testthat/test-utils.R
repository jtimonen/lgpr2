test_that("extend_df2 works", {
  dat <- testdata
  t <- c(1, 3)
  w <- 40
  df <- extend_df2(dat, t, "t", w, "weight", "id")
  expect_true(all(df$w == 40))
  num_ids <- length(unique(dat$id))
  expect_true(nrow(df) == num_ids * length(t))
})

test_that("extend_df2 works", {
  dat <- testdata
  t <- c(1, 3)
  w <- 40
  df <- extend_df2(dat, t, "t", w, "weight", "id")
  expect_true(all(df$w == 40))
  num_ids <- length(unique(dat$id))
  expect_true(nrow(df) == num_ids * length(t))
})
