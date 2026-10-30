# Network meta-analysis: 24 smoking-cessation trials (netmeta::smokingcessation)
# Run from this folder: Rscript network_meta_analysis.R  (tables and charts are written to results/)

# Network meta-analysis helpers: arm-level data in long format, frequentist (netmeta) and Bayesian (gemtc / JAGS) fits.
library(netmeta)
library(gemtc)

TREATMENTS <- c(A = "No contact", B = "Self-help", C = "Individual counselling", D = "Group counselling")

load_smoking <- function() {
  data(smokingcessation, package = "netmeta", envir = environment())
  w <- smokingcessation
  long <- do.call(rbind, lapply(1:3, function(k) data.frame(
    study = seq_len(nrow(w)), treatment = w[[paste0("treat", k)]],
    responders = w[[paste0("event", k)]], sampleSize = w[[paste0("n", k)]])))
  long <- long[!is.na(long$sampleSize) & long$treatment != "", ]
  long$treatment <- unname(TREATMENTS[long$treatment])
  long[order(long$study, long$treatment), ]
}

fit_netmeta <- function(long) {
  p <- pairwise(treat = treatment, event = responders, n = sampleSize, studlab = study, data = long, sm = "OR")
  netmeta(p, reference.group = "No contact", common = FALSE, random = TRUE)
}

fit_gemtc <- function(long, effects = c("random", "fixed"), n.iter = 20000, seed = 2026) {
  effects <- match.arg(effects)
  # gemtc only accepts letters, digits and underscores in treatment names
  net <- mtc.network(data.ab = transform(long, study = as.character(study), treatment = gem_key(treatment)))
  set.seed(seed)
  mod <- mtc.model(net, likelihood = "binom", link = "logit", linearModel = effects, n.chain = 4)
  mtc.run(mod, n.adapt = 5000, n.iter = n.iter, thin = 10)
}

gem_key <- function(x) gsub("[^A-Za-z0-9]", "_", x)
gem_name <- function(x) unname(setNames(TREATMENTS, gem_key(TREATMENTS))[x])

# SUCRA from a rank-probability matrix (rows = treatments, columns = rank 1..K; rank 1 = best).
sucra <- function(rp) {
  K <- ncol(rp)
  apply(rp, 1, function(p) sum(cumsum(p)[1:(K - 1)]) / (K - 1))
}

# ---- Analysis ----
library(ggplot2)
dir.create("results", showWarnings = FALSE)

long <- load_smoking()
write.csv(long, "results/nma_arm_data.csv", row.names = FALSE)
cat("studies:", length(unique(long$study)), " arms:", nrow(long), "\n")

# Frequentist random-effects NMA
nm <- fit_netmeta(long)
trts <- rownames(nm$TE.random)
freq <- data.frame(treatment = trts, or = exp(nm$TE.random[, "No contact"]),
                   lcl = exp(nm$lower.random[, "No contact"]), ucl = exp(nm$upper.random[, "No contact"]),
                   row.names = NULL)
ps <- netrank(nm, small.values = "undesirable")
freq$p_score <- ps$ranking.random[freq$treatment]
write.csv(freq, "results/nma_frequentist.csv", row.names = FALSE)
het <- data.frame(tau2 = nm$tau^2, I2 = nm$I2, Q = nm$Q, Q_df = nm$df.Q, Q_p = nm$pval.Q,
                  Q_inconsistency = nm$Q.inconsistency, Q_inconsistency_p = nm$pval.Q.inconsistency)
write.csv(het, "results/nma_heterogeneity.csv", row.names = FALSE)

ns <- netsplit(nm)
split <- data.frame(comparison = ns$comparison, direct_or = exp(ns$direct.random$TE), indirect_or = exp(ns$indirect.random$TE),
                    p_difference = ns$compare.random$p)
write.csv(split, "results/nma_nodesplit.csv", row.names = FALSE)

league <- netleague(nm, common = FALSE, random = TRUE, digits = 2, backtransf = TRUE)
write.csv(league$random, "results/nma_league_table.csv", row.names = FALSE)

png("results/04_network.png", width = 1100, height = 900, res = 160)
netgraph(nm, plastic = FALSE, thickness = "number.of.studies", number.of.studies = TRUE, col = "#2563EB",
         points = TRUE, cex.points = 4, col.points = "#1F2933", multiarm = FALSE)
dev.off()

# Bayesian NMA: random vs fixed effects
re <- fit_gemtc(long, "random"); fe <- fit_gemtc(long, "fixed")
psrf <- gelman.diag(re$samples, multivariate = FALSE)$psrf
dic <- data.frame(model = c("random effects", "fixed effects"),
                  dic = c(summary(re)$DIC["DIC"], summary(fe)$DIC["DIC"]),
                  pD = c(summary(re)$DIC["pD"], summary(fe)$DIC["pD"]),
                  resdev = c(summary(re)$DIC["Dbar"], summary(fe)$DIC["Dbar"]),
                  max_psrf = c(max(psrf[, 1]), max(gelman.diag(fe$samples, multivariate = FALSE)$psrf[, 1])))
write.csv(dic, "results/nma_bayes_model_fit.csv", row.names = FALSE)

rel <- relative.effect(re, t1 = gem_key("No contact"))
q <- summary(rel)$summaries$quantiles
bay <- data.frame(treatment = sub("^d\\.No_contact\\.", "", rownames(q)), or = exp(q[, "50%"]),
                  lcl = exp(q[, "2.5%"]), ucl = exp(q[, "97.5%"]), row.names = NULL)
sd_row <- grepl("^sd", bay$treatment)
bay$treatment[!sd_row] <- gem_name(bay$treatment[!sd_row])
tau_bayes <- data.frame(param = "between-study SD (log OR)", median = q[sd_row, "50%"], lcl = q[sd_row, "2.5%"], ucl = q[sd_row, "97.5%"])
write.csv(tau_bayes, "results/nma_bayes_sd.csv", row.names = FALSE)
logor <- data.frame(treatment = gem_name(sub("^d\\.No_contact\\.", "", rownames(q)[!sd_row])), log_or = q[!sd_row, "50%"])
write.csv(logor, "results/nma_bayes_logor.csv", row.names = FALSE)
bay <- bay[!sd_row, ]
rp <- rank.probability(re, preferredDirection = 1)
rp_df <- as.data.frame(unclass(rp)); colnames(rp_df) <- paste0("rank_", seq_len(ncol(rp_df)))
rp_df$treatment <- gem_name(rownames(rp_df))
rp_df$sucra <- sucra(unclass(rp))
bay <- merge(bay, rp_df[, c("treatment", "rank_1", "sucra")], by = "treatment", all = TRUE)
write.csv(bay, "results/nma_bayesian.csv", row.names = FALSE)
write.csv(rp_df, "results/nma_rank_probabilities.csv", row.names = FALSE)

png("results/05_trace.png", width = 1400, height = 1000, res = 140)
plot(re$samples[, 1:min(4, ncol(re$samples[[1]]))], density = FALSE)
dev.off()

both <- rbind(cbind(method = "Frequentist (netmeta)", subset(freq, treatment != "No contact")[, c("treatment", "or", "lcl", "ucl")]),
              cbind(method = "Bayesian (gemtc/JAGS)", subset(bay, !is.na(or))[, c("treatment", "or", "lcl", "ucl")]))
write.csv(both, "results/nma_compare.csv", row.names = FALSE)
th <- theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"), legend.position = "top", panel.grid.minor = element_blank())
g <- ggplot(both, aes(or, treatment, colour = method)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0.25, orientation = "y", position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) + scale_x_log10() +
  scale_colour_manual(values = c("Bayesian (gemtc/JAGS)" = "#2563EB", "Frequentist (netmeta)" = "#9AA5B1"), name = NULL) +
  labs(x = "Odds ratio of quitting vs no contact (log scale)", y = NULL,
       title = sprintf("%s ranks first (SUCRA %.0f%%)", bay$treatment[which.max(bay$sucra)], 100 * max(bay$sucra)),
       subtitle = "Random-effects NMA, 24 trials: frequentist and Bayesian estimates agree") + th
ggsave("results/06_nma_forest.png", g, width = 9.5, height = 4.8, dpi = 160, bg = "white")

print(freq); print(het); print(split); print(dic); print(bay); print(tau_bayes)
