# Survival helpers: Kaplan-Meier summaries, Cox models, parametric extrapolation and mean survival.
library(survival)
library(flexsurv)

DISTS <- c(exponential = "exp", weibull = "weibull", gompertz = "gompertz",
           lognormal = "lnorm", loglogistic = "llogis", gengamma = "gengamma")

load_colon <- function() {
  d <- subset(survival::colon, etype == 2)            # death records only
  d$years <- d$time / 365.25
  d$age65 <- factor(ifelse(d$age >= 65, ">= 65", "< 65"), levels = c("< 65", ">= 65"))
  d$sexf <- factor(d$sex, 0:1, c("Female", "Male"))
  d$nodes4 <- factor(d$node4, 0:1, c("<= 4 nodes", "> 4 nodes"))
  d
}

km_summary <- function(d) {
  fit <- survfit(Surv(years, status) ~ rx, data = d)
  tab <- summary(fit)$table
  s5 <- summary(fit, times = 5)
  rm5 <- summary(fit, rmean = 5)$table[, "rmean"]
  data.frame(arm = sub("rx=", "", rownames(tab)), n = tab[, "records"], deaths = tab[, "events"],
             median_years = tab[, "median"], median_lcl = tab[, "0.95LCL"], median_ucl = tab[, "0.95UCL"],
             surv_5y = s5$surv, rmst_5y = unname(rm5), row.names = NULL)
}

hr_table <- function(fit, label) {
  s <- summary(fit)$conf.int
  data.frame(model = label, term = rownames(s), hr = s[, 1], lcl = s[, 3], ucl = s[, 4],
             p = summary(fit)$coefficients[, "Pr(>|z|)"], row.names = NULL)
}

# Fit every candidate distribution to one arm; return fits plus an AIC/BIC table.
fit_parametric <- function(d_arm) {
  fits <- lapply(DISTS, function(dist) flexsurvreg(Surv(years, status) ~ 1, data = d_arm, dist = dist))
  ic <- data.frame(model = names(fits), aic = sapply(fits, AIC), bic = sapply(fits, BIC), row.names = NULL)
  ic$delta_aic <- ic$aic - min(ic$aic)
  list(fits = fits, ic = ic[order(ic$aic), ])
}

surv_curve <- function(fit, times) {
  summary(fit, t = times, type = "survival", ci = FALSE, tidy = TRUE)$est
}

# Mean survival = area under S(t) to the horizon (trapezoid on a fine grid).
mean_survival <- function(fit, horizon = 30, step = 0.01) {
  t <- seq(0, horizon, by = step)
  s <- surv_curve(fit, t)
  sum(diff(t) * (head(s, -1) + tail(s, -1)) / 2)
}

subgroup_hr <- function(d, var, arm = "Lev+5FU", ref = "Obs") {
  d <- droplevels(subset(d, rx %in% c(ref, arm)))
  d$rx <- relevel(d$rx, ref)
  out <- lapply(levels(d[[var]]), function(lv) {
    s <- summary(coxph(Surv(years, status) ~ rx, data = d[d[[var]] == lv, ]))$conf.int
    data.frame(subgroup = var, level = lv, n = sum(d[[var]] == lv), hr = s[1, 1], lcl = s[1, 3], ucl = s[1, 4])
  })
  inter <- anova(coxph(as.formula(paste("Surv(years, status) ~ rx +", var)), data = d),
                 coxph(as.formula(paste("Surv(years, status) ~ rx *", var)), data = d))
  res <- do.call(rbind, out)
  res$p_interaction <- inter[["Pr(>|Chi|)"]][2]
  res
}
