source(file.path("..", "..", "R", "survival.R"))

test_that("colon death records load with 929 patients in 3 arms", {
  d <- load_colon()
  expect_equal(nrow(d), 929)
  expect_equal(as.vector(table(d$rx)), c(315, 310, 304))
  expect_true(all(d$years >= 0))
})

test_that("mean survival matches the closed form for an exponential model", {
  d <- load_colon()
  fit <- flexsurvreg(Surv(years, status) ~ 1, data = subset(d, rx == "Obs"), dist = "exp")
  rate <- fit$res["rate", "est"]
  expect_equal(mean_survival(fit, horizon = 30), (1 - exp(-rate * 30)) / rate, tolerance = 1e-3)
})

test_that("KM summary covers every patient and subgroup HRs come with an interaction test", {
  d <- load_colon()
  expect_equal(sum(km_summary(d)$n), 929)
  sg <- subgroup_hr(d, "sexf")
  expect_equal(nrow(sg), 2)
  expect_true(all(sg$lcl < sg$hr & sg$hr < sg$ucl))
  expect_true(sg$p_interaction[1] > 0 && sg$p_interaction[1] < 1)
})
