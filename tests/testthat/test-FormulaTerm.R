test_that("stancode creation for a GPTerm works", {
  x <- GPTerm$new("x", NA)
  expect_output(print(x), "f_gp_x")
  code1 <- x$stancode_tdata(datanames = c("OBS"))
  code2 <- x$stancode_tdata(datanames = c("OBS", "PRD"))
  expect_gt(nchar(code1), 100) # 198?
  expect_gt(nchar(code2), nchar(code1))
})
