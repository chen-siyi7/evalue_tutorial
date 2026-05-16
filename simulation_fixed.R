# =============================================================================
# Fixed-sample simulation for the e-value tutorial.
#
# Reproduces Table 3: empirical rejection rates at alpha = 0.05 for K = 4
# correlated one-sided z-tests under three dependence structures
# (independence, equicorrelation rho = 0.5, two-block correlation rho = 0.5)
# and three mean vectors (null; sparse alternative (3,0,0,0); dense (1,1,1,1)).
#
# Methods:
#   - Fisher's combination (assumes independence)
#   - Stouffer's combination, assuming independence
#   - Stouffer's combination, using the true Sigma
#   - Bonferroni
#   - Cauchy combination test (ACAT)
#   - Mean e-value at kappa = 0.5, uniform weights
#   - Mean e-value, outer mixture over kappa in {0.1, 0.3, 0.5, 0.7, 0.9}
#
# M = 10^5 replicates per scenario; MC standard error at p = 0.05 is ~0.0007.
#
# Requires: mvtnorm
# =============================================================================

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("Please install the 'mvtnorm' package: install.packages('mvtnorm').")
}

set.seed(202605)

K     <- 4
alpha <- 0.05
M     <- 1e5

make_sigma <- function(kind, rho = 0.5, K = 4) {
  if (kind == "indep") {
    return(diag(K))
  }
  if (kind == "equicorr") {
    return((1 - rho) * diag(K) + rho * matrix(1, K, K))
  }
  if (kind == "block") {
    S <- diag(K)
    S[1, 2] <- S[2, 1] <- rho
    S[3, 4] <- S[4, 3] <- rho
    return(S)
  }
  stop(sprintf("unknown dependence kind: %s", kind))
}

simulate_z <- function(mu, Sigma, M) {
  # rows are replicates, columns are streams
  mvtnorm::rmvnorm(n = M, mean = mu, sigma = Sigma)
}

fisher_reject <- function(P, alpha, K) {
  stat  <- -2 * rowSums(log(P))
  pcomb <- pchisq(stat, df = 2 * K, lower.tail = FALSE)
  pcomb < alpha
}

stouffer_indep_reject <- function(Z, alpha, K) {
  Zstat <- rowSums(Z) / sqrt(K)
  pcomb <- pnorm(Zstat, lower.tail = FALSE)
  pcomb < alpha
}

stouffer_correct_reject <- function(Z, Sigma, alpha) {
  Zstat <- rowSums(Z) / sqrt(sum(Sigma))
  pcomb <- pnorm(Zstat, lower.tail = FALSE)
  pcomb < alpha
}

bonferroni_reject <- function(P, alpha, K) {
  pcomb <- pmin(K * apply(P, 1, min), 1)
  pcomb < alpha
}

cct_reject <- function(P, alpha) {
  # Liu and Xie (2020) Cauchy combination, equal weights;
  # the null distribution of the combined statistic is approximately Cauchy.
  eps <- 1e-15
  Pc  <- pmin(pmax(P, eps), 1 - eps)
  Tstat <- rowMeans(tan((0.5 - Pc) * pi))
  pcomb <- 0.5 - atan(Tstat) / pi
  pcomb < alpha
}

mean_evalue_at_kappa <- function(P, kappa) {
  # f_kappa(u) = kappa * u^(kappa-1); valid calibrator with uniform weights
  rowMeans(kappa * P^(kappa - 1))
}

mean_evalue_reject <- function(P, alpha, kappa = 0.5) {
  mean_evalue_at_kappa(P, kappa) >= 1 / alpha
}

outer_mix_reject <- function(P, alpha, kappa_grid = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  Ebars <- sapply(kappa_grid, function(k) mean_evalue_at_kappa(P, k))
  Eouter <- rowMeans(Ebars)
  Eouter >= 1 / alpha
}

run_scenario <- function(name, kind, mu, M, alpha, K, rho = 0.5) {
  Sigma <- make_sigma(kind, rho = rho, K = K)
  Z     <- simulate_z(mu, Sigma, M)
  P     <- pnorm(Z, lower.tail = FALSE)

  rates <- c(
    "Fisher"          = mean(fisher_reject(P, alpha, K)),
    "Stouffer-IND"    = mean(stouffer_indep_reject(Z, alpha, K)),
    "Stouffer-COR"    = mean(stouffer_correct_reject(Z, Sigma, alpha)),
    "Bonferroni"      = mean(bonferroni_reject(P, alpha, K)),
    "Cauchy"          = mean(cct_reject(P, alpha)),
    "MeanE-0.5"       = mean(mean_evalue_reject(P, alpha, kappa = 0.5)),
    "MeanE-outer"     = mean(outer_mix_reject(P, alpha))
  )

  cat(sprintf("\n=== %s: dependence=%s, rho=%g ===\n", name, kind, rho))
  cat(sprintf("  mu = (%s)\n", paste(format(mu, nsmall = 2), collapse = ", ")))
  for (nm in names(rates)) {
    r  <- rates[[nm]]
    se <- sqrt(r * (1 - r) / M)
    cat(sprintf("    %-14s: %.4f  (MC s.e. %.4f)\n", nm, r, se))
  }
  invisible(rates)
}

# ----- scenarios reproducing Table 3 -----------------------------------------

zero   <- rep(0, K)
sparse <- c(3, 0, 0, 0)
dense  <- c(1, 1, 1, 1)

results <- list(
  null_indep   = run_scenario("NULL (indep)",        "indep",    zero,   M, alpha, K),
  null_eq05    = run_scenario("NULL (equicorr 0.5)", "equicorr", zero,   M, alpha, K, rho = 0.5),
  null_block   = run_scenario("NULL (block 0.5)",    "block",    zero,   M, alpha, K, rho = 0.5),
  sparse_indep = run_scenario("SPARSE (indep)",      "indep",    sparse, M, alpha, K),
  sparse_eq05  = run_scenario("SPARSE (eq 0.5)",     "equicorr", sparse, M, alpha, K, rho = 0.5),
  dense_indep  = run_scenario("DENSE (indep)",       "indep",    dense,  M, alpha, K),
  dense_eq05   = run_scenario("DENSE (eq 0.5)",      "equicorr", dense,  M, alpha, K, rho = 0.5)
)

cat(sprintf(
  "\n\nSUMMARY\nM = %d, K = %d, alpha = %g\nMC s.e. at p=0.05 is ~%.4f\n",
  M, K, alpha, sqrt(alpha * (1 - alpha) / M)
))
