# Statistical analysis plan

Written before the analysis was run. Deviations are listed at the end.

## 1. Objectives
**A. Survival (single trial, individual patient data).** For patients with resected stage C colon cancer, estimate the
effect of adjuvant levamisole + 5-fluorouracil (Lev+5FU) and levamisole alone (Lev) on overall survival against
observation (Obs), and extrapolate survival beyond trial follow-up to estimate the mean life-years gained, as an
economic model submitted to NICE would need.

**B. Network meta-analysis (aggregate data).** Compare four smoking-cessation interventions (no contact, self-help,
individual counselling, group counselling) using all direct and indirect trial evidence, and rank them.

## 2. Data
- A: `survival::colon` (Moertel et al., 1990; 929 patients, 3 arms; death records only, `etype == 2`). Time in days, converted to years.
- B: `netmeta::smokingcessation` (Hasselblad 1998; 24 studies, 4 treatments; arm-level counts of quitters).

## 3. Analyses: survival
1. Kaplan–Meier estimates by arm; median survival with 95% CI; log-rank test (3 arms and Lev+5FU vs Obs).
2. Cox proportional hazards: unadjusted, then adjusted for pre-specified prognostic factors (age, sex, obstruction,
   perforation, adherence, extent, more than 4 positive nodes, differentiation). Report HR (95% CI).
3. Proportional hazards check: Schoenfeld residuals (`cox.zph`), global and per covariate, and log(−log) plot.
4. Parametric extrapolation, fitted separately per arm (NICE DSU TSD 14): exponential, Weibull, Gompertz, log-normal,
   log-logistic, generalised gamma. Select by AIC/BIC **and** visual fit **and** clinical plausibility of the long-term tail.
5. Restricted mean survival time (RMST) at 5 years from KM, and mean survival over 30 years from the selected models;
   life-years gained = difference between arms.
6. Subgroups (pre-specified, exploratory): sex, age (< 65 / ≥ 65), more than 4 positive nodes. HR for Lev+5FU vs Obs with
   interaction p-values. No multiplicity adjustment; interpreted as hypothesis-generating only.

## 4. Analyses: network meta-analysis
1. Network plot; check connectivity.
2. Frequentist random-effects NMA on the log odds ratio scale (`netmeta`), reference = no contact. Report ORs (95% CI),
   league table, P-scores, between-study heterogeneity (tau², I²).
3. Inconsistency: design-by-treatment interaction test and node-splitting (direct vs indirect).
4. Bayesian NMA (`gemtc`, JAGS = BUGS language): binomial likelihood, logit link, fixed and random effects; vague priors;
   4 chains, 20,000 iterations after 5,000 burn-in. Convergence: Gelman–Rubin PSRF < 1.05 and trace plots.
   Model choice by DIC. Report ORs (95% credible intervals), rank probabilities and SUCRA.
5. Compare frequentist and Bayesian estimates.

## 5. Software and quality
R 4.5; packages `survival`, `flexsurv`, `netmeta`, `gemtc`, `rjags`. Seeds fixed. Functions unit-tested with `testthat`
against known values; tests run in GitHub Actions. All numbers in the README come from `outputs/`.

## 6. Deviations from plan
1. **Adjusted Cox model is complete-case** (906 of 929 patients): 23 patients have missing tumour differentiation. Not imputed.
2. **Proportional hazards:** holds for treatment (p = 0.28) but not globally (p < 0.001), driven by differentiation,
   obstruction and > 4 nodes. The treatment HR is reported as planned; a stratified or time-varying sensitivity analysis
   for those covariates was not run and is listed as a next step.
3. **No general-population mortality cap** on the extrapolations (needs ONS life tables by age and sex). Planned for a later version.
4. R 4.5 was used instead of R 4.4 (the version available for Apple Silicon at the time).
