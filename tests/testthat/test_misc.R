test_that("Check for when no arguments are passed to Trace_plot", {
  expect_error(Trace_plot(), "variable and chain must be provided\n")
})

test_that("Title & labels for Trace_plot are as expected", {
  LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                Cens = cancer[, 2], X = cancer[, 3:11])

  p <- Trace_plot(1, LN)
  expect_equal(as.character(p$labels[1]),
               "Trace Plot for Variable beta.Intercept")
  expect_equal(as.character(p$labels[2]), "Iteration")
  expect_equal(as.character(p$labels[3]), "Value")
})

test_that("BASSLINE_convert returns expected value", {
  Time <- c(5, 15, 15)
  Cens <- c(1, 0, 1)
  experiment <- as.factor(c("chem1", "chem2", "chem3"))
  age <- c(15, 35, 20)

  df <- data.frame(Time, Cens, experiment, age)
  converted <- BASSLINE_convert(df)
  expected <- cbind(c(5, 15, 15), c(1, 0, 1), c(1, 0, 0), c(0, 1, 0),
                    c(0, 0, 1), c(15, 35, 20))
  colnames(expected) <- c("Time", "Cens", "experiment.chem1",
                          "experiment.chem2", "experiment.chem3", "age")
  row.names(expected) <- seq(3)
  expect_equal(converted, as.matrix(expected))
})
