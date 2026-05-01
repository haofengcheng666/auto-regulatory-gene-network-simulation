# Post-processing for the adaptive-tau tolerance pilot.
# Reads the per-method RDS files written by scripts/pilot_sweep.R and
# produces:
#   outputs/pilot_summary_N<N>.csv    -- one row per setting
#   outputs/pilot_traj_mean_N<N>.csv  -- ensemble mean trajectories
#   outputs/pilot_traj_sd_N<N>.csv    -- ensemble sd trajectories
#   outputs/pilot_table_N<N>.tex      -- LaTeX table for Appendix A.5
#
# RMSE rows include 95% paired bootstrap CIs (B resamples of replicate
# indices, shared across method and SSA so the matched-seed structure is
# preserved).
#
# Usage:
#   Rscript scripts/pilot_postprocess.R                       # N=100, B=1000
#   Rscript scripts/pilot_postprocess.R --N=1000 --B=1000     # full size
#   Rscript scripts/pilot_postprocess.R --N=1000 outputs_test

args <- commandArgs(trailingOnly = TRUE)
N_pilot <- 100L
B_boot  <- 1000L
output_dir <- "outputs"
for (a in args) {
  if (grepl("^--N=", a)) {
    N_pilot <- as.integer(sub("^--N=", "", a))
  } else if (grepl("^--B=", a)) {
    B_boot <- as.integer(sub("^--B=", "", a))
  } else if (nzchar(a) && !startsWith(a, "--")) {
    output_dir <- a
  }
}
cat(sprintf("pilot_postprocess.R  N = %d  B = %d  output_dir = %s\n",
            N_pilot, B_boot, output_dir))

rrc_grid <- c(0.05, 0.10, 0.15, 0.20)

ssa_path <- file.path(output_dir, sprintf("pilot_partial_ssa_N%d.rds", N_pilot))
stopifnot(file.exists(ssa_path))
ssa <- readRDS(ssa_path)

adapt <- list()
for (rrc in rrc_grid) {
  p <- file.path(output_dir,
                 sprintf("pilot_partial_rrc_%.2f_N%d.rds", rrc, N_pilot))
  stopifnot(file.exists(p))
  adapt[[sprintf("%.2f", rrc)]] <- readRDS(p)
}

t_grid  <- ssa$t_grid
n_grid  <- length(t_grid)
P2_ssa  <- ssa$P2
rt_ssa  <- ssa$runtime_s
mean_ssa <- rowMeans(P2_ssa, na.rm = TRUE)
sd_ssa   <- apply(P2_ssa, 1, sd, na.rm = TRUE)

rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))

paired_bootstrap_ci <- function(P2_method, P2_ref, B, master_seed = 42L) {
  set.seed(master_seed)
  vals_mean <- numeric(B)
  vals_sd   <- numeric(B)
  N <- ncol(P2_method)
  for (b in seq_len(B)) {
    idx <- sample.int(N, N, replace = TRUE)
    sub_a <- P2_method[, idx, drop = FALSE]
    sub_r <- P2_ref[, idx, drop = FALSE]
    m_a <- rowMeans(sub_a, na.rm = TRUE)
    s_a <- apply(sub_a, 1, sd, na.rm = TRUE)
    m_r <- rowMeans(sub_r, na.rm = TRUE)
    s_r <- apply(sub_r, 1, sd, na.rm = TRUE)
    vals_mean[b] <- rmse(m_a, m_r)
    vals_sd[b]   <- rmse(s_a, s_r)
    if (b %% max(1L, B %/% 20L) == 0L || b == B) {
      cat(sprintf("    bootstrap %d/%d\n", b, B)); flush.console()
    }
  }
  list(rmse_mean_lo = unname(quantile(vals_mean, 0.025, na.rm = TRUE)),
       rmse_mean_hi = unname(quantile(vals_mean, 0.975, na.rm = TRUE)),
       rmse_sd_lo   = unname(quantile(vals_sd,   0.025, na.rm = TRUE)),
       rmse_sd_hi   = unname(quantile(vals_sd,   0.975, na.rm = TRUE)))
}

rows <- list()
mean_traj <- data.frame(time = t_grid, SSA = mean_ssa)
sd_traj   <- data.frame(time = t_grid, SSA = sd_ssa)

for (rrc in rrc_grid) {
  cat(sprintf("---- post-process rrc = %.2f ----\n", rrc))
  blk <- adapt[[sprintf("%.2f", rrc)]]
  P2_a <- blk$P2
  rt_a <- blk$runtime_s
  ne_a <- blk$n_events
  m_a  <- rowMeans(P2_a, na.rm = TRUE)
  s_a  <- apply(P2_a, 1, sd, na.rm = TRUE)
  mean_traj[[sprintf("Adaptive_rrc%.2f", rrc)]] <- m_a
  sd_traj[[sprintf("Adaptive_rrc%.2f",   rrc)]] <- s_a

  diffs <- rt_a - rt_ssa
  paired_se <- sd(diffs) / sqrt(length(diffs))
  bs <- paired_bootstrap_ci(P2_a, P2_ssa, B = B_boot)

  rows[[sprintf("rrc%.2f", rrc)]] <- data.frame(
    rrc = rrc, n_pilot = N_pilot,
    mean_runtime_s = mean(rt_a), sd_runtime_s = sd(rt_a),
    paired_diff_vs_ssa_s = mean(diffs),
    paired_diff_lo = mean(diffs) - 1.96 * paired_se,
    paired_diff_hi = mean(diffs) + 1.96 * paired_se,
    rmse_mean_P2 = rmse(m_a, mean_ssa),
    rmse_mean_lo = bs$rmse_mean_lo,
    rmse_mean_hi = bs$rmse_mean_hi,
    rmse_sd_P2 = rmse(s_a, sd_ssa),
    rmse_sd_lo = bs$rmse_sd_lo,
    rmse_sd_hi = bs$rmse_sd_hi,
    mean_n_events = mean(ne_a, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

ssa_summary <- data.frame(
  rrc = NA_real_, n_pilot = N_pilot,
  mean_runtime_s = mean(rt_ssa), sd_runtime_s = sd(rt_ssa),
  paired_diff_vs_ssa_s = 0, paired_diff_lo = 0, paired_diff_hi = 0,
  rmse_mean_P2 = 0, rmse_mean_lo = 0, rmse_mean_hi = 0,
  rmse_sd_P2 = 0, rmse_sd_lo = 0, rmse_sd_hi = 0,
  mean_n_events = NA_real_,
  stringsAsFactors = FALSE)

summary_df <- do.call(rbind, c(list(ssa_summary), rows))

out_csv <- file.path(output_dir, sprintf("pilot_summary_N%d.csv", N_pilot))
write.csv(summary_df, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")
write.csv(mean_traj,
          file.path(output_dir, sprintf("pilot_traj_mean_N%d.csv", N_pilot)),
          row.names = FALSE)
write.csv(sd_traj,
          file.path(output_dir, sprintf("pilot_traj_sd_N%d.csv", N_pilot)),
          row.names = FALSE)

# LaTeX table
fmt_d <- function(x, k = 2) formatC(x, digits = k, format = "f")
fmt_signed <- function(x, k = 2) {
  s <- ifelse(x >= 0, "+", "")
  paste0(s, fmt_d(x, k))
}
tex <- c(
  "% Auto-generated by scripts/pilot_postprocess.R; do not hand-edit.",
  sprintf("%% N_pilot = %d, B_boot = %d, smfsb v1.5, adaptivetau v2.3.2",
          N_pilot, B_boot),
  "\\begin{table}[ht]",
  "\\centering",
  sprintf("\\caption[Adaptive-tolerance pilot]{Adaptive $\\tau$-leaping tolerance pilot at $N_{\\mathrm{pilot}} = %d$ replicates per setting, with paired SSA reference under the same per-replicate seed vector. The paired runtime difference is reported as $\\overline{\\Delta} = T^{(\\mathrm{ada})} - T^{(\\mathrm{SSA})}$ with a $95\\%%$ normal-approximation CI; the RMSE summaries are computed against the SSA reference on the common grid $T = \\{0,1,\\dots,500\\}$, with $95\\%%$ paired bootstrap CIs ($B = %d$ resamples of replicate indices, shared between adaptive and SSA so the matched-seed structure is preserved).}",
          N_pilot, B_boot),
  "\\label{tab:pilot}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4.0pt}%",
  "\\begin{tabular}{@{}cccrcccc@{}}",
  "\\toprule",
  "\\texttt{rrc} & Mean runtime (s) & sd (s) & $\\overline{\\Delta}$ vs SSA (s) & $95\\%$ CI & $\\mathrm{RMSE}_{\\mathrm{mean}}(P_2)$ & $\\mathrm{RMSE}_{\\mathrm{sd}}(P_2)$ & mean $n_{\\text{ev}}$ \\\\",
  "\\midrule"
)
ssa_row <- summary_df[is.na(summary_df$rrc), ]
tex <- c(tex,
  sprintf("SSA (ref) & %s & %s & %s & --- & %s & %s & --- \\\\",
          fmt_d(ssa_row$mean_runtime_s),
          fmt_d(ssa_row$sd_runtime_s),
          fmt_d(0), fmt_d(0), fmt_d(0)))
for (rrc in rrc_grid) {
  r <- summary_df[!is.na(summary_df$rrc) & summary_df$rrc == rrc, ]
  ev_str <- if (is.na(r$mean_n_events)) "---" else
    paste0(formatC(r$mean_n_events / 1000, digits = 1, format = "f"), "k")
  tex <- c(tex,
    sprintf("$%.2f$ & %s & %s & %s & $[%s,\\, %s]$ & $%s\\ [%s,\\, %s]$ & $%s\\ [%s,\\, %s]$ & %s \\\\",
            rrc,
            fmt_d(r$mean_runtime_s),
            fmt_d(r$sd_runtime_s),
            fmt_signed(r$paired_diff_vs_ssa_s),
            fmt_d(r$paired_diff_lo),
            fmt_d(r$paired_diff_hi),
            fmt_d(r$rmse_mean_P2),
            fmt_d(r$rmse_mean_lo),
            fmt_d(r$rmse_mean_hi),
            fmt_d(r$rmse_sd_P2),
            fmt_d(r$rmse_sd_lo),
            fmt_d(r$rmse_sd_hi),
            ev_str))
}
tex <- c(tex, "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}")
out_tex <- file.path(output_dir, sprintf("pilot_table_N%d.tex", N_pilot))
writeLines(tex, out_tex)
cat("Wrote", out_tex, "\n")

cat("\n----- summary -----\n")
print(summary_df, row.names = FALSE, digits = 4)
