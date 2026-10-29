source(file.path("..", "..", "R", "nma.R"))

test_that("smoking network has 24 studies, 50 arms and 4 treatments", {
  long <- load_smoking()
  expect_equal(length(unique(long$study)), 24)
  expect_equal(nrow(long), 50)
  expect_setequal(unique(long$treatment), unname(TREATMENTS))
})

test_that("gemtc-safe names round-trip", {
  expect_equal(gem_name(gem_key(TREATMENTS)), unname(TREATMENTS))
})

test_that("SUCRA is 1 for a certain best and 0 for a certain worst", {
  rp <- diag(3); rownames(rp) <- c("best", "middle", "worst")
  expect_equal(unname(sucra(rp)), c(1, 0.5, 0))
})

test_that("frequentist NMA puts every active treatment above no contact", {
  nm <- fit_netmeta(load_smoking())
  te <- nm$TE.random[, "No contact"]
  expect_true(all(te[names(te) != "No contact"] > 0))
})

test_that("Bayesian results reproduce the NICE DSU TSD2 smoking example", {
  f <- file.path("..", "..", "outputs", "nma_bayes_logor.csv")
  skip_if_not(file.exists(f))
  got <- setNames(read.csv(f)$log_or, read.csv(f)$treatment)
  expect_equal(unname(got[c("Self-help", "Individual counselling", "Group counselling")]), c(0.49, 0.84, 1.10), tolerance = 0.05)
  sd <- read.csv(file.path("..", "..", "outputs", "nma_bayes_sd.csv"))$median
  expect_equal(sd, 0.82, tolerance = 0.05)
})
