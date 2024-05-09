test_that("creating a Stan model and Stan code for empty model works", {
  sm <- StanModel$new(compile = FALSE)
  sc <- sm$create_stancode()
  expect_true(inherits(sm, "StanModel"))
  expect_gt(nchar(sc), 100)
})
