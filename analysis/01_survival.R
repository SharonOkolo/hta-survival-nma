# Objective A: survival analysis and extrapolation (see SAP.md section 3). Run from the repo root.
source("R/survival.R")
library(ggplot2)
set.seed(2026)
dir.create("outputs", showWarnings = FALSE); dir.create("screenshots", showWarnings = FALSE)

d <- load_colon()

# 1. Kaplan-Meier and log-rank
km <- km_summary(d)
write.csv(km, "outputs/km_summary.csv", row.names = FALSE)
lr3 <- survdiff(Surv(years, status) ~ rx, data = d)
d2 <- droplevels(subset(d, rx %in% c("Obs", "Lev+5FU")))
lr2 <- survdiff(Surv(years, status) ~ rx, data = d2)
logrank <- data.frame(comparison = c("3 arms", "Lev+5FU vs Obs"),
                      chisq = c(lr3$chisq, lr2$chisq),
                      df = c(2, 1),
                      p = c(pchisq(lr3$chisq, 2, lower.tail = FALSE), pchisq(lr2$chisq, 1, lower.tail = FALSE)))
write.csv(logrank, "outputs/logrank.csv", row.names = FALSE)

kmfit <- survfit(Surv(years, status) ~ rx, data = d)
km_pts <- data.frame(arm = rep(sub("rx=", "", names(kmfit$strata)), kmfit$strata), time = kmfit$time,
                     surv = kmfit$surv, lower = kmfit$lower, upper = kmfit$upper)
write.csv(km_pts, "outputs/km_points.csv", row.names = FALSE)

# 2. Cox models
cox_u <- coxph(Surv(years, status) ~ rx, data = d)
covs <- "age + sex + obstruct + perfor + adhere + factor(extent) + node4 + factor(differ)"
cox_a <- coxph(as.formula(paste("Surv(years, status) ~ rx +", covs)), data = d)
cox_tab <- rbind(hr_table(cox_u, "unadjusted"), hr_table(cox_a, "adjusted"))
write.csv(cox_tab, "outputs/cox.csv", row.names = FALSE)
n_complete <- cox_a$n

# 3. Proportional hazards check
zph <- cox.zph(cox_a)
ph <- data.frame(term = rownames(zph$table), chisq = zph$table[, "chisq"], p = zph$table[, "p"], row.names = NULL)
write.csv(ph, "outputs/ph_test.csv", row.names = FALSE)

# 4. Parametric extrapolation per arm
arms <- c("Obs", "Lev+5FU")
pf <- lapply(arms, function(a) fit_parametric(subset(d, rx == a)))
names(pf) <- arms
ic <- do.call(rbind, lapply(arms, function(a) cbind(arm = a, pf[[a]]$ic)))
write.csv(ic, "outputs/parametric_fit.csv", row.names = FALSE)

times <- seq(0, 30, by = 0.25)
curves <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(names(DISTS), function(m)
  data.frame(arm = a, model = m, time = times, surv = surv_curve(pf[[a]]$fits[[m]], times))))))
write.csv(curves, "outputs/extrapolation_curves.csv", row.names = FALSE)

# 5. Mean survival (30-year horizon) for every model; the base case is chosen in the README with reasons.
ms <- do.call(rbind, lapply(names(DISTS), function(m) data.frame(
  model = m,
  obs = mean_survival(pf[["Obs"]]$fits[[m]]),
  lev5fu = mean_survival(pf[["Lev+5FU"]]$fits[[m]]))))
ms$life_years_gained <- ms$lev5fu - ms$obs
ms$surv_20y_obs <- sapply(names(DISTS), function(m) surv_curve(pf[["Obs"]]$fits[[m]], 20))
ms$surv_20y_lev5fu <- sapply(names(DISTS), function(m) surv_curve(pf[["Lev+5FU"]]$fits[[m]], 20))
write.csv(ms, "outputs/mean_survival.csv", row.names = FALSE)

# Parameters for the Excel model (Weibull and log-normal, flexsurv's parameterisation)
params <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(c("weibull", "lognormal", "loglogistic"), function(m) {
  r <- pf[[a]]$fits[[m]]$res
  data.frame(arm = a, model = m, param = rownames(r), est = r[, "est"], row.names = NULL)
}))))
write.csv(params, "outputs/parametric_params.csv", row.names = FALSE)

# 6. Subgroups
sg <- do.call(rbind, lapply(c("sexf", "age65", "nodes4"), function(v) subgroup_hr(d, v)))
sg$subgroup <- c(sexf = "Sex", age65 = "Age", nodes4 = "Positive nodes")[sg$subgroup]
write.csv(sg, "outputs/subgroups.csv", row.names = FALSE)

# Figures
pal <- c("Obs" = "#9AA5B1", "Lev" = "#F0B429", "Lev+5FU" = "#2563EB")
th <- theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"), legend.position = "top", panel.grid.minor = element_blank())
g1 <- ggplot(km_pts, aes(time, surv, colour = arm)) + geom_step(linewidth = 0.9) +
  scale_colour_manual(values = pal, name = NULL) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Years since randomisation", y = "Overall survival",
       title = sprintf("Lev+5FU cut the risk of death by %.0f%% vs observation", 100 * (1 - exp(coef(cox_u))["rxLev+5FU"])),
       subtitle = "Kaplan-Meier, colon cancer trial (n = 929)") + th
ggsave("screenshots/01_km.png", g1, width = 9, height = 5.2, dpi = 160, bg = "white")

base <- subset(curves, model %in% c("weibull", "lognormal", "loglogistic", "gompertz"))
g2 <- ggplot() +
  geom_line(data = base, aes(time, surv, colour = model), linewidth = 0.8) +
  geom_step(data = subset(km_pts, arm %in% arms), aes(time, surv), colour = "black", linewidth = 0.6) +
  facet_wrap(~arm) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Years", y = "Overall survival", colour = "Model",
       title = "Same data, very different futures: the tail drives life-years gained",
       subtitle = "Kaplan-Meier (black) vs parametric extrapolations to 30 years") + th
ggsave("screenshots/02_extrapolation.png", g2, width = 10, height = 5.2, dpi = 160, bg = "white")

sg$label <- paste(sg$subgroup, sg$level, sep = ": ")
g3 <- ggplot(sg, aes(hr, reorder(label, nrow(sg):1))) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0.2, orientation = "y", colour = "#2563EB") + geom_point(size = 3, colour = "#2563EB") +
  scale_x_log10() + labs(x = "Hazard ratio, Lev+5FU vs Obs (log scale)", y = NULL,
  title = sprintf("Benefit looks larger in men (interaction p = %.2f)", sg$p_interaction[sg$subgroup == "Sex"][1]),
  subtitle = "Exploratory, unadjusted for multiplicity: hypothesis-generating only") + th
ggsave("screenshots/03_subgroups.png", g3, width = 9, height = 4.5, dpi = 160, bg = "white")

cat("complete cases in adjusted Cox:", n_complete, "\n")
print(km); print(logrank); print(subset(cox_tab, grepl("^rx", term))); print(ph); print(ic); print(ms); print(sg)
