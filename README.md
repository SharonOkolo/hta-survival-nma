# Survival extrapolation and network meta-analysis

Two analyses behind most NICE health technology assessments, done in R on public trial data. Read
**[Project_Report.docx](Project_Report.docx)** for the full write-up: the problem, methods, results, charts and the
statistical analysis plan.

| File | What it is |
|---|---|
| `Project_Report.docx` | The project explained, with results and charts |
| `survival_analysis.R` | Kaplan–Meier, Cox models, proportional hazards check, 6 parametric extrapolations, subgroups |
| `network_meta_analysis.R` | Frequentist (netmeta) and Bayesian (gemtc + JAGS) network meta-analysis |
| `survival_model.xlsx` | Excel survival model built from the R estimates: pick a model, see life-years gained |

**Headline results**
- Levamisole + 5-FU cut the risk of death by 31% in a 929-patient colon cancer trial (HR 0.69, 95% CI 0.55–0.87).
- Extrapolated to 30 years: 4.6 life-years gained (log-normal base case), but 2.9 to 5.8 depending on the model.
- Across 24 smoking-cessation trials, group counselling ranks first (OR about 3.0 vs no contact, 95% CrI about 1.3–7.3),
  matching the published NICE DSU worked example.

**Run it** (R with JAGS installed; packages: survival, flexsurv, netmeta, gemtc, rjags, ggplot2)
```bash
Rscript survival_analysis.R
Rscript network_meta_analysis.R
```
