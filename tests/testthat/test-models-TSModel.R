test_that("generating only the Stan code works", {
  code1 <- stancode_ts(y ~ gp(x), print = FALSE)
  expect_match(code1, "Term f_gp_x")
  expect_gt(nchar(code1), 20)
})


test_that("creating a TSModel works", {
  m <- TSModel$new(hello ~ gp(foo) + gp(foo, bar))
  expect_output(m$print(), "TSModel")
  expect_output(cat(m$term_list$latex()), "alpha")
  expect_output(m$term_list$terms[[1]]$print())
  expect_equal(m$y_var, "hello")
  expect_equal(m$term_list$length(), 3)
  tn <- m$term_list$terms[[1]]$name_long()
  expect_true(is.character(tn))
  expect_gt(nchar(tn), 3)
})

test_that("creating the Stan data works", {
  m <- TSModel$new(hello ~ gp(foo) + gp(foo, bar), compile = F, id_var = "bar")
  a <- data.frame(
    foo = c(-100, 2, 3, 4),
    hello = c(0, 3, 2, 1),
    bar = as.factor(c(1, 1, 2, 2))
  )
  expect_error(m$fit(), "you need to call compile")
  sbf <- 2.7
  confs <- m$term_list$fill_term_confs(num_bf = 2, scale_bf = sbf)
  m$term_list$set_transforms(a)
  sd <- m$create_standata(a, term_confs = confs)
  sd <- sd$stan_data
  expect_equal(sd$n_LON, 4)
  expect_equal(sd$B_foo, 2)
  expect_equal(max(sd$dat_foo_unit_LON), 1)
  expect_equal(min(sd$dat_foo_unit_LON), -1)
  expect_equal(sd$L_foo, sbf)
  expect_equal(m$id_var, "bar")
})


test_that("fitting a model and plotting function draws work", {
  m <- TSModel$new(hello ~ gp(foo) + gp(foo, bar), delta = 0.05)
  a <- data.frame(
    foo = c(1, 2, 3, 4),
    hello = c(0, 3, 2, 1),
    bar = as.factor(c(1, 1, 2, 2))
  )
  num_bf <- 8
  fit <- m$fit(data = a, refresh = 0, num_bf = num_bf)
  B_foo <- fit$term_confs[["f_gp_foo"]]$num_bf
  expect_equal(B_foo, num_bf)
  expect_equal(m$get_delta(), 0.05)
  fd <- fit$function_draws()
  f1 <- fit$function_draws("f_gp_foo")
  f2 <- fit$function_draws("f_gp_fooXbar")
  f3 <- f1 + f2
  f4 <- f3 - f1
  expect_output(fd$print(), "FunctionDraws")
  expect_output(f1$print(), "FunctionDraws")
  expect_output(f2$print(), "FunctionDraws")
  expect_output(f3$print(), "FunctionDraws")
  expect_output(f4$print(), "FunctionDraws")
  expect_s3_class(fd$plot(), "ggplot")
  expect_s3_class(f1$plot(), "ggplot")
  expect_s3_class(f2$plot(), "ggplot")
  expect_s3_class(f3$plot(), "ggplot")
  expect_s3_class(f4$plot(), "ggplot")
  expect_s3_class(fit$plot(), "ggplot")
  plt <- fit$plot(f_reference = rep(1, 4))
  ll <- fit$loglik()
  ep <- fit$measurement_error()
  expect_equal(length(ll), nrow(a))
  expect_equal(length(ep), nrow(a))
  expect_s3_class(plt, "ggplot")

  df <- fit$fit_quality_summary()
  expect_equal(nrow(df), 1)
  expect_equal(ncol(df), 2) # colnames c("id", "error_perc")
})

test_that("a gp example works", {
  r <- example(
    formula = "y ~ gp(x)",
    iter_warmup = 500, iter_sampling = 500, chains = 1
  )
  p1 <- r$plot()
  p2 <- (r$function_draws(data_scale = FALSE) - r$function_draws("f_gp_x"))$plot()
  p3 <- (r$function_draws("f_gp_x") + r$function_draws("f_baseline_id"))$plot()
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  p4 <- r$predict()$plot(predictive = FALSE, capped = FALSE)
  expect_s3_class(p4, "ggplot")
})

test_that("simplest model with empty formula (only grouped offset) works", {
  a <- TSModel$new(y ~ .)
  r <- a$fit(testdata)
  expect_equal(ncol(r$function_draws()$get_input()), 1)
})
