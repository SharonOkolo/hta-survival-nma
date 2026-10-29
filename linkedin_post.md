# LinkedIn post draft

Six survival models. All of them fit the trial data well. Over 30 years, they disagree by almost 3 life-years.

That gap is why statisticians matter in health technology assessment, so I built the two analyses behind most NICE
submissions, end to end in R:

🩺 The problem: trials stop after a few years and rarely compare every treatment, yet decisions are made on lifetime
benefit across all the options.

⚙️ The process: I wrote the statistical analysis plan first, then:
• Survival analysis of a 929-patient colon cancer trial: Kaplan–Meier, Cox models, a proportional hazards check and
  6 parametric extrapolations per arm
• Frequentist and Bayesian network meta-analysis of 24 smoking-cessation trials (netmeta, and gemtc with JAGS)
• An Excel model generated from the R estimates, checked to match R within 1%
• 10 automated tests, including one that reproduces the published NICE DSU worked example

📊 The results:
• Levamisole + 5-FU cut the risk of death by 31% (HR 0.69, 95% CI 0.55–0.87)
• Base case: 4.6 life-years gained, but 2.9 to 5.8 depending on the model, which is why scenario analyses matter
• Group counselling ranked first for quitting smoking (OR 2.99 vs no contact; 71% probability of being best)
• A subgroup "finding" (bigger benefit in men, p = 0.04) that I reported as hypothesis-generating, not a reason to treat differently

🔗 Code: https://github.com/SharonOkolo/hta-survival-nma
🔗 App: [Streamlit link]

#Biostatistics #HealthEconomics #HTA #RStats #SurvivalAnalysis #MetaAnalysis #NICE

---
**Attach (in this order):** `screenshots/02_extrapolation.png`, `screenshots/01_km.png`, `screenshots/06_nma_forest.png`,
`screenshots/07_app_survival.png`.
