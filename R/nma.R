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
