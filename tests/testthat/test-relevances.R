test_that("relevance computation works", {
  a <- example(chains = 1)
  b <- a$relevances()
  expect_equal(length(b), 3)
})
