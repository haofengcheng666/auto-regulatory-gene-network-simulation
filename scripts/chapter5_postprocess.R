# Across-replicate paired bootstrap on the chapter 5 RMSE summaries.
# Reads the five partial_*.rds files written by benchmark_N1000_final.R
# and writes outputs/chapter5_summary_with_bootstrap_N<N>.csv with the
# matched-seed paired bootstrap CIs that go into Table 5.2.
#
# Same bootstrap as scripts/pilot_postprocess.R: at each iteration draw
# one set of N replicate indices and apply it to BOTH the method matrix
# and the SSA reference matrix, so the matched-seed pairing carries
# through; recompute the bootstrapped ensemble mean / sd trajectories
# and the RMSE summaries; take percentile CIs over B iterations.
#
# Usage:
#   Rscript scripts/chapter5_postprocess.R               # outputs/, B=1000
#   Rscript scripts/chapter5_postprocess.R --B=2000 outputs_test

args <- commandArgs(trailingOnly = TRUE)
B_boot      <- 1000L
master_seed <- 42L
output_dir  <- "outputs"
for (a in args) {
  if (grepl("^--B=", a)) {
    B_boot <- as.integer(sub("^--B=", "", a))
  } else if (grepl("^--seed=", a)) {
    master_seed <- as.integer(sub("^--seed=", "", a))
  } else if (nzchar(a) && !startsWith(a, "--")) {
    output_dir <- a
  }
}
cat(sprintf("chapter5_postprocess.R  B = %d  seed = %d  output_dir = %s\n",
            B_boot, master_seed, output_dir))

required <- c("partial_ssa.rds", "partial_adaptive.rds",
              "partial_fixed_0.25.rds", "partial_fixed_0.5.rds",
              "partial_fixed_1.0.rds")
missing <- required[!file.exists(file.path(output_dir, required))]
if (length(missing) > 0L) {
  stop("Missing partials in ", output_dir, ": ",
       paste(missing, collapse = ", "),
       ". Run benchmark_N1000_final.R first.")
}

partials <- lapply(required, function(f) readRDS(file.path(output_dir, f)))
names(partials) <- gsub("^partial_|\\.rds$", "", required)

t_grid <- partials[[1]]$t_grid
N_rep  <- ncol(partials[[1]]$P2_matrix)

P2_ssa   <- partials$ssa$P2_matrix
mean_ssa <- rowMeans(P2_ssa, na.rm = TRUE)
sd_ssa   <- apply(P2_ssa, 1, sd, na.rm = TRUE)

rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))

paired_bootstrap_ci <- function(P2_method, P2_ref, B, seed) {
  set.seed(seed)
  vals_mean <- numeric(B)
  vals_sd   <- numeric(B)
  N <- ncol(P2_method)
  for (b in seq_len(B)) {
    idx <- sample.int(N, N, replace = TRUE)
    sub_a <- P2_method[, idx, drop = FALSE]
    sub_r <- P2_ref[,    idx, drop = FALSE]
    m_a <- rowMeans(sub_a, na.rm = TRUE)
    s_a <- apply(sub_a, 1, sd, na.rm = TRUE)
    m_r <- rowMeans(sub_r, na.rm = TRUE)
    s_r <- apply(sub_r, 1, sd, na.rm = TRUE)
    vals_mean[b] <- rmse(m_a, m_r)
    vals_sd[b]   <- rmse(s_a, s_r)
    if (b %% max(1L, B %/% 10L) == 0L || b == B) {
      cat(sprintf("    bootstrap %d/%d\n", b, B)); flush.console()
    }
  }
  list(rmse_mean_lo = unname(quantile(vals_mean, 0.025, na.rm = TRUE)),
       rmse_mean_hi = unname(quantile(vals_mean, 0.975, na.rm = TRUE)),
       rmse_sd_lo   = unname(quantile(vals_sd,   0.025, na.rm = TRUE)),
       rmse_sd_hi   = unname(quantile(vals_sd,   0.975, na.rm = TRUE)))
}

method_label <- function(k) {
  switch(k,
         ssa        = "SSA",
         adaptive   = "Adaptive tau (rrc=0.15)",
         fixed_0.25 = "Fixed tau=0.25",
         fixed_0.5  = "Fixed tau=0.5",
         fixed_1.0  = "Fixed tau=1.0")
}
tau_label <- function(k) {
  switch(k, ssa = "-", adaptive = "0.15",
         fixed_0.25 = "0.25", fixed_0.5 = "0.5", fixed_1.0 = "1")
}

rows <- list()
# SSA reference: zero RMSE by construction, no bootstrap needed
rows[["ssa"]] <- data.frame(
  method = method_label("ssa"), tau_or_rrc = tau_label("ssa"), n_rep = N_rep,
  rmse_mean_P2 = 0, rmse_mean_lo = 0, rmse_mean_hi = 0,
  rmse_sd_P2   = 0, rmse_sd_lo   = 0, rmse_sd_hi   = 0,
  bootstrap_B  = B_boot, bootstrap_seed = master_seed,
  stringsAsFactors = FALSE
)

for (k in c("adaptive", "fixed_0.25", "fixed_0.5", "fixed_1.0")) {
  cat(sprintf("---- bootstrap %s ----\n", k))
  P2_m <- partials[[k]]$P2_matrix
  m    <- rowMeans(P2_m, na.rm = TRUE)
  s    <- apply(P2_m, 1, sd, na.rm = TRUE)
  bs   <- paired_bootstrap_ci(P2_m, P2_ssa,
                              B = B_boot, seed = master_seed)
  rows[[k]] <- data.frame(
    method = method_label(k), tau_or_rrc = tau_label(k), n_rep = N_rep,
    rmse_mean_P2 = rmse(m, mean_ssa),
    rmse_mean_lo = bs$rmse_mean_lo,
    rmse_mean_hi = bs$rmse_mean_hi,
    rmse_sd_P2   = rmse(s, sd_ssa),
    rmse_sd_lo   = bs$rmse_sd_lo,
    rmse_sd_hi   = bs$rmse_sd_hi,
    bootstrap_B  = B_boot, bootstrap_seed = master_seed,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, rows)
out_csv <- file.path(output_dir,
                     sprintf("chapter5_summary_with_bootstrap_N%d.csv", N_rep))
write.csv(summary_df, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")

cat("\n----- summary -----\n")
print(summary_df, row.names = FALSE, digits = 4)
