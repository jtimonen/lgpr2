test_that("stanmodel_functions works", {
  sm <- stanmodel_functions("gp/basisfun")
  expect_true(inherits(sm, "CmdStanModel"))
})
