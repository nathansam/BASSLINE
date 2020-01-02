test_that("Check for when no arguments are passed to Trace_Plot",{
  expect_error(Trace_Plot(), "obs and chain must be provided\n")
})

test_that("Title & labels for Trace_Plot are as expected",{
  LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                Cens = cancer[,2], X = cancer[,3:11])

  p <- Trace_Plot(1, LN)
  expect_equal(as.character(p$labels[1]), "Trace Plot for Observation 1")
  expect_equal(as.character(p$labels[2]), "Iteration")
  expect_equal(as.character(p$labels[3]), "Value")
})
