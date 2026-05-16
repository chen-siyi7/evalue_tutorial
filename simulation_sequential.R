# =============================================================================
# Sequential simulation for the e-value tutorial.
#
# Reproduces Table 4: empirical rejection rates under repeated optional
# inspection. K = 4 correlated streams accumulate data in T = 10 equally
# spaced batches of n_batch = 20 observations per stream. The analyst tests
# at every interim and rejects at the first crossing.
#
# Methods (all at nominal alpha = 0.05):
#   - Fisher's combination, applied at every look without alpha-spending
#   - Stouffer's combination with the true Sigma, at every look
#   - Mean e-process built from one likelihood-ratio test martingale per
#     stream against the point alternative mu_alt = 0.20; rejection on first
#     crossing of 1/alpha = 20.
#
# Cross-stream dependence is equicorrelated with rho = 0.5 (or independence,
# for the reference column).
#
# M = 5e4 replicates per scenario.
#
# Requires: mvtnorm
# =============================================================================

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("Please install the 'mvtnorm' package: install.packages('mvtnorm').")
}

set.seed(202606)

K           <- 4
T_looks     <- 10
n_batch     <- 20
alpha       <- 0.05
M           <- 5e4
mu_alt_e    <- 0.20    # the LR-test martingale's betting point

make_corr <- function(rho, K = 4) {
  (1 - rho) * diag(K) + rho * matrix(1, K, K)
}

simulate_streams <- function(mu_true, Sigma, M, T_looks, n_batch) {
  # Returns:
  #   Z         : (M, T_looks, K) cumulative z-statistics through look t
  #   cumS_hist : (M, T_looks, K) cumulative sum of observations per stream
  L <- t(chol(Sigma))                # lower-triangular Cholesky factor
  cumS      <- matrix(0, M, K)
  Z         <- array(0, dim = c(M, T_looks, K))
  cumS_hist <- array(0, dim = c(M, T_looks, K))

  for (t in seq_len(T_looks)) {
    # draw a batch of n_batch obs per replicate per stream
    batch_total <- matrix(0, M, K)
    for (i in seq_len(n_batch)) {
      white  <- matrix(rnorm(M * K), M, K)
      sample <- white %*% t(L)                            # (M, K)
      sample <- sweep(sample, 2, mu_true, "+")
      batch_total <- batch_total + sample
    }
    cumS <- cumS + batch_total
    n_t  <- t * n_batch
    Z[,         t, ] <- cumS / sqrt(n_t)
    cumS_hist[, t, ] <- cumS
  }
  list(Z = Z, cumS_hist = cumS_hist)
}

fisher_rejects <- function(Z, alpha, K) {
  # Z: (M, T, K). Test at every look without alpha-spending; reject at
  # first crossing.
  dims <- dim(Z)
  M_   <- dims[1]; T_ <- dims[2]
  P    <- pmin(pmax(pnorm(Z, lower.tail = FALSE), 1e-300), 1)
  # collapse the third axis (streams) to a per-look chi-square statistic
  stat  <- -2 * apply(log(P), c(1, 2), sum)              # (M, T)
  pcomb <- pchisq(stat, df = 2 * K, lower.tail = FALSE)  # (M, T)
  apply(pcomb < alpha, 1, any)
}

stouffer_correct_rejects <- function(Z, Sigma, alpha) {
  # Stouffer with the true Sigma at every look.
  var_sum <- sum(Sigma)
  Zsum    <- apply(Z, c(1, 2), sum)                     # (M, T)
  Zstat   <- Zsum / sqrt(var_sum)
  pcomb   <- pnorm(Zstat, lower.tail = FALSE)
  apply(pcomb < alpha, 1, any)
}

mean_eprocess_rejects <- function(cumS_hist, n_batch, mu_alt, alpha) {
  # E_{k,t} = exp(mu_alt * S_{k,t} - n_t * mu_alt^2 / 2). Under H_0 each is
  # a nonneg martingale starting at 1; the mean is an e-process by Theorem 2.
  dims <- dim(cumS_hist)
  T_   <- dims[2]
  K_   <- dims[3]

  # n_t along the time axis
  n_vec <- n_batch * seq_len(T_)
  # log-component (M, T, K): subtract (n_t * mu_alt^2 / 2) along axis 2
  logE <- mu_alt * cumS_hist
  logE <- sweep(logE, 2, n_vec * (mu_alt^2) / 2, "-")
  E    <- exp(logE)
  Ebar <- apply(E, c(1, 2), mean)                       # (M, T)
  apply(Ebar >= 1 / alpha, 1, any)
}

run_scenario <- function(name, mu_vec, rho, M, T_looks, n_batch, K,
                         alpha, mu_alt_e) {
  Sigma <- make_corr(rho, K = K)
  sims  <- simulate_streams(as.numeric(mu_vec), Sigma, M, T_looks, n_batch)

  fish  <- mean(fisher_rejects(sims$Z, alpha, K))
  stou  <- mean(stouffer_correct_rejects(sims$Z, Sigma, alpha))
  eproc <- mean(mean_eprocess_rejects(sims$cumS_hist, n_batch, mu_alt_e, alpha))

  cat(sprintf("\n=== %s ===\n", name))
  cat(sprintf("  rho=%g, K=%d, T=%d looks, n_batch=%d, M=%d\n",
              rho, K, T_looks, n_batch, M))
  cat(sprintf("  mu = (%s)\n", paste(format(mu_vec), collapse = ", ")))
  for (entry in list(c("Fisher (any look)",              fish),
                     c("Stouffer-true Sigma (any look)", stou),
                     c("Mean e-process (any look)",      eproc))) {
    nm <- entry[[1]]
    r  <- as.numeric(entry[[2]])
    se <- sqrt(r * (1 - r) / M)
    cat(sprintf("    %-32s: %.4f  (MC s.e. %.4f)\n", nm, r, se))
  }
}

# ----- scenarios reproducing Table 4 -----------------------------------------

zero       <- c(0, 0, 0, 0)
sparse_seq <- c(0.15, 0, 0, 0)        # one stream with weak per-obs drift
dense_seq  <- c(0.07, 0.07, 0.07, 0.07)

run_scenario("SEQUENTIAL NULL (rho=0.5)",   zero,        0.5,
             M, T_looks, n_batch, K, alpha, mu_alt_e)
run_scenario("SEQUENTIAL SPARSE (rho=0.5)", sparse_seq,  0.5,
             M, T_looks, n_batch, K, alpha, mu_alt_e)
run_scenario("SEQUENTIAL DENSE  (rho=0.5)", dense_seq,   0.5,
             M, T_looks, n_batch, K, alpha, mu_alt_e)
# Reference: independence under the null
run_scenario("SEQUENTIAL NULL (indep)",     zero,        0.0,
             M, T_looks, n_batch, K, alpha, mu_alt_e)
