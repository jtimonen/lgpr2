test_that("relevance computation works", {
  a <- example2(chains = 1)
  b <- a$relevances()
  expect_equal(length(b), 3)
})
