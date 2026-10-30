# Survival analysis and extrapolation: colon cancer trial (survival::colon)
# Run from this folder: Rscript survival_analysis.R  (tables and charts are written to results/)

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

# ---- Analysis ----
library(ggplot2)
set.seed(2026)
dir.create("results", showWarnings = FALSE)

d <- load_colon()

# 1. Kaplan-Meier and log-rank
km <- km_summary(d)
write.csv(km, "results/km_summary.csv", row.names = FALSE)
lr3 <- survdiff(Surv(years, status) ~ rx, data = d)
d2 <- droplevels(subset(d, rx %in% c("Obs", "Lev+5FU")))
lr2 <- survdiff(Surv(years, status) ~ rx, data = d2)
logrank <- data.frame(comparison = c("3 arms", "Lev+5FU vs Obs"),
                      chisq = c(lr3$chisq, lr2$chisq),
                      df = c(2, 1),
                      p = c(pchisq(lr3$chisq, 2, lower.tail = FALSE), pchisq(lr2$chisq, 1, lower.tail = FALSE)))
write.csv(logrank, "results/logrank.csv", row.names = FALSE)

kmfit <- survfit(Surv(years, status) ~ rx, data = d)
km_pts <- data.frame(arm = rep(sub("rx=", "", names(kmfit$strata)), kmfit$strata), time = kmfit$time,
                     surv = kmfit$surv, lower = kmfit$lower, upper = kmfit$upper)
write.csv(km_pts, "results/km_points.csv", row.names = FALSE)

# 2. Cox models
cox_u <- coxph(Surv(years, status) ~ rx, data = d)
covs <- "age + sex + obstruct + perfor + adhere + factor(extent) + node4 + factor(differ)"
cox_a <- coxph(as.formula(paste("Surv(years, status) ~ rx +", covs)), data = d)
cox_tab <- rbind(hr_table(cox_u, "unadjusted"), hr_table(cox_a, "adjusted"))
write.csv(cox_tab, "results/cox.csv", row.names = FALSE)
n_complete <- cox_a$n

# 3. Proportional hazards check
zph <- cox.zph(cox_a)
ph <- data.frame(term = rownames(zph$table), chisq = zph$table[, "chisq"], p = zph$table[, "p"], row.names = NULL)
write.csv(ph, "results/ph_test.csv", row.names = FALSE)

# 4. Parametric extrapolation per arm
arms <- c("Obs", "Lev+5FU")
pf <- lapply(arms, function(a) fit_parametric(subset(d, rx == a)))
names(pf) <- arms
ic <- do.call(rbind, lapply(arms, function(a) cbind(arm = a, pf[[a]]$ic)))
write.csv(ic, "results/parametric_fit.csv", row.names = FALSE)

times <- seq(0, 30, by = 0.25)
curves <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(names(DISTS), function(m)
  data.frame(arm = a, model = m, time = times, surv = surv_curve(pf[[a]]$fits[[m]], times))))))
write.csv(curves, "results/extrapolation_curves.csv", row.names = FALSE)

# 5. Mean survival (30-year horizon) for every model; the base case is chosen in the README with reasons.
ms <- do.call(rbind, lapply(names(DISTS), function(m) data.frame(
  model = m,
  obs = mean_survival(pf[["Obs"]]$fits[[m]]),
  lev5fu = mean_survival(pf[["Lev+5FU"]]$fits[[m]]))))
ms$life_years_gained <- ms$lev5fu - ms$obs
ms$surv_20y_obs <- sapply(names(DISTS), function(m) surv_curve(pf[["Obs"]]$fits[[m]], 20))
ms$surv_20y_lev5fu <- sapply(names(DISTS), function(m) surv_curve(pf[["Lev+5FU"]]$fits[[m]], 20))
write.csv(ms, "results/mean_survival.csv", row.names = FALSE)

# Parameters for the Excel model (Weibull and log-normal, flexsurv's parameterisation)
params <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(c("weibull", "lognormal", "loglogistic"), function(m) {
  r <- pf[[a]]$fits[[m]]$res
  data.frame(arm = a, model = m, param = rownames(r), est = r[, "est"], row.names = NULL)
}))))
write.csv(params, "results/parametric_params.csv", row.names = FALSE)

# 6. Subgroups
sg <- do.call(rbind, lapply(c("sexf", "age65", "nodes4"), function(v) subgroup_hr(d, v)))
sg$subgroup <- c(sexf = "Sex", age65 = "Age", nodes4 = "Positive nodes")[sg$subgroup]
write.csv(sg, "results/subgroups.csv", row.names = FALSE)

# Figures
pal <- c("Obs" = "#9AA5B1", "Lev" = "#F0B429", "Lev+5FU" = "#2563EB")
th <- theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"), legend.position = "top", panel.grid.minor = element_blank())
g1 <- ggplot(km_pts, aes(time, surv, colour = arm)) + geom_step(linewidth = 0.9) +
  scale_colour_manual(values = pal, name = NULL) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Years since randomisation", y = "Overall survival",
       title = sprintf("Lev+5FU cut the risk of death by %.0f%% vs observation", 100 * (1 - exp(coef(cox_u))["rxLev+5FU"])),
       subtitle = "Kaplan-Meier, colon cancer trial (n = 929)") + th
ggsave("results/01_km.png", g1, width = 9, height = 5.2, dpi = 160, bg = "white")

base <- subset(curves, model %in% c("weibull", "lognormal", "loglogistic", "gompertz"))
g2 <- ggplot() +
  geom_line(data = base, aes(time, surv, colour = model), linewidth = 0.8) +
  geom_step(data = subset(km_pts, arm %in% arms), aes(time, surv), colour = "black", linewidth = 0.6) +
  facet_wrap(~arm) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Years", y = "Overall survival", colour = "Model",
       title = "Same data, very different futures: the tail drives life-years gained",
       subtitle = "Kaplan-Meier (black) vs parametric extrapolations to 30 years") + th
ggsave("results/02_extrapolation.png", g2, width = 10, height = 5.2, dpi = 160, bg = "white")

sg$label <- paste(sg$subgroup, sg$level, sep = ": ")
g3 <- ggplot(sg, aes(hr, reorder(label, nrow(sg):1))) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0.2, orientation = "y", colour = "#2563EB") + geom_point(size = 3, colour = "#2563EB") +
  scale_x_log10() + labs(x = "Hazard ratio, Lev+5FU vs Obs (log scale)", y = NULL,
  title = sprintf("Benefit looks larger in men (interaction p = %.2f)", sg$p_interaction[sg$subgroup == "Sex"][1]),
  subtitle = "Exploratory, unadjusted for multiplicity: hypothesis-generating only") + th
ggsave("results/03_subgroups.png", g3, width = 9, height = 4.5, dpi = 160, bg = "white")

cat("complete cases in adjusted Cox:", n_complete, "\n")
print(km); print(logrank); print(subset(cox_tab, grepl("^rx", term))); print(ph); print(ic); print(ms); print(sg)
