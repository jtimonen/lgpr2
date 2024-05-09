test_that("stanmodel_functions works", {
  sm <- stanmodel_functions("hazard/integrate")
  expect_true(inherits(sm, "CmdStanModel"))
})
