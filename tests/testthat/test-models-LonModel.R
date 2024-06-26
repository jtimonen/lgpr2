test_that("creating a LonModel works", {
  m <- LonModel$new(hello ~ gp(foo) + gp(foo, bar) + offset(id))
  expect_output(m$print(), "LonModel")
  expect_output(m$term_list$terms[[1]]$print())
  expect_equal(m$y_var, "hello")
  expect_equal(m$term_list$length(), 3)
  tn <- m$term_list$terms[[1]]$name_long()
  expect_true(is.character(tn))
  expect_gt(nchar(tn), 3)
})

test_that("creating the Stan data works", {
  m <- LonModel$new(hello ~ gp(foo) + gp(foo, bar) + offset(bar), compile = F)
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
})


test_that("fitting a model and plotting function draws work", {
  m <- LonModel$new(hello ~ gp(foo) + gp(foo, bar))
  a <- data.frame(
    foo = c(1, 2, 3, 4),
    hello = c(0, 3, 2, 1),
    bar = as.factor(c(1, 1, 2, 2))
  )
  num_bf <- 8
  fit <- m$fit(data = a, refresh = 0, num_bf = num_bf)
  B_foo <- fit$term_confs[["f_gp_foo"]]$num_bf
  expect_equal(B_foo, num_bf)
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
  expect_equal(length(ll), nrow(a))
  expect_s3_class(plt, "ggplot")
})

test_that("a gp example works with predict", {
  r <- example(
    formula = "y ~ gp(x)",
    iter_warmup = 500, iter_sampling = 500, chains = 1
  )
  d <- r$diagnose()
  expect_equal(length(d), 4)
  red <- r$reduce()
  expect_equal(red$selected[1], "p_noise")
  p1 <- r$plot()
  p2 <- (r$function_draws() - r$function_draws("f_gp_x"))$plot()
  p3 <- (r$function_draws() + r$function_draws("f_gp_x"))$plot()
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  p4 <- r$predict()$plot(predictive = FALSE)
  expect_s3_class(p4, "ggplot")

  # Predict works
  p <- r$predict_time(t_test = seq(0, 12, by = 0.1), t_var = "x")
  expect_s3_class(p$plot(), "ggplot")

  # Should not work outside [-L, L]
  expect_error(
    {
      p <- r$predict_time(t_test = seq(0, 22, by = 0.2), t_var = "x")
    },
    " GP approximation not valid for this input"
  )
})
