# Make the chapter 5 figures from the N=1000 benchmark CSVs.
# Run from the benchmark/ folder. No re-simulation here, only read CSV.
#   Rscript scripts/make_chapter5_figures.R                       # outputs/ -> figures/
#   Rscript scripts/make_chapter5_figures.R outputs_test figures_test

data_dir <- "outputs"
out_dir  <- "figures"
.cli <- commandArgs(trailingOnly = TRUE)
if (length(.cli) >= 1 && nzchar(.cli[1])) data_dir <- .cli[1]
if (length(.cli) >= 2 && nzchar(.cli[2])) out_dir  <- .cli[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# simple base R colours, one per method
cols <- c("black", "blue", "orange", "red", "darkgreen")

method_keys   <- c("SSA", "Adaptive", "Fixed_tau_0.25",
                   "Fixed_tau_0.5", "Fixed_tau_1")
method_labels <- c("SSA", "Adaptive tau",
                   "Fixed tau=0.25", "Fixed tau=0.5", "Fixed tau=1.0")

summary_file <- file.path(data_dir, "benchmark_summary_N1000.csv")
if (!file.exists(summary_file)) {
  stop("Missing ", summary_file, ", run the benchmark first.")
}
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)

# build a key column so we can match the summary rows to method_keys
row_key <- function(r) {
  m   <- trimws(r[["method"]])
  tau <- trimws(as.character(r[["tau"]]))
  if (m == "SSA") return("SSA")
  if (startsWith(m, "Adaptive")) return("Adaptive")
  if (startsWith(m, "Fixed")) {
    n <- as.numeric(tau)
    return(paste0("Fixed_tau_", format(n, trim = TRUE, nsmall = 0)))
  }
  m
}
summary_df$key <- vapply(seq_len(nrow(summary_df)),
                         function(i) row_key(summary_df[i, ]),
                         character(1))

plot_runtime <- function() {
  idx  <- match(method_keys, summary_df$key)
  vals <- summary_df$mean_runtime_pipeline_s[idx]
  sds  <- summary_df$sd_runtime_pipeline_s[idx]

  pdf(file.path(out_dir, "fig_runtime.pdf"), width = 6, height = 3.5)
  par(mar = c(4, 8, 1, 2), las = 1)
  rev_vals <- rev(vals); rev_sds <- rev(sds)
  bp <- barplot(rev_vals, names.arg = rev(method_labels), horiz = TRUE,
                col = rev(cols), xlab = "Mean runtime per replicate (s)",
                xlim = c(0, max(vals + sds) * 1.15))
  # +/- 1 sd whiskers on each bar
  arrows(x0 = rev_vals - rev_sds, y0 = bp,
         x1 = rev_vals + rev_sds, y1 = bp,
         code = 3, angle = 90, length = 0.04, lwd = 1, col = "grey20")
  dev.off()
}

plot_safeguard <- function() {
  keys <- c("Fixed_tau_0.25", "Fixed_tau_0.5", "Fixed_tau_1")
  idx  <- match(keys, summary_df$key)
  vals <- summary_df$mean_safeguard_triggers[idx]
  labs <- method_labels[3:5]
  cls  <- cols[3:5]

  pdf(file.path(out_dir, "fig_safeguard.pdf"), width = 6, height = 2.8)
  par(mar = c(4, 8, 1, 2), las = 1)
  barplot(rev(vals), names.arg = rev(labs), horiz = TRUE,
          col = rev(cls),
          xlab = "Mean unsafe nominal proposals per replicate",
          xlim = c(0, max(vals) * 1.15))
  dev.off()
}

# one helper for both mean and sd trajectory plots
trajectory_plot <- function(df, ylab, out_file,
                            legend_corner = "ul_inset") {
  t <- df$time
  pdf(out_file, width = 6.5, height = 4)
  par(mar = c(4, 4, 1, 1))
  yl <- c(0, max(unlist(df[, -1]), na.rm = TRUE) * 1.05)
  plot(t, df$SSA, type = "l", col = cols[1], lwd = 2,
       xlab = "Time", ylab = ylab, ylim = yl)
  if ("Adaptive"       %in% names(df)) lines(t, df$Adaptive,       col = cols[2], lwd = 1.5)
  if ("Fixed_tau_0.25" %in% names(df)) lines(t, df$Fixed_tau_0.25, col = cols[3], lwd = 1)
  if ("Fixed_tau_0.5"  %in% names(df)) lines(t, df$Fixed_tau_0.5,  col = cols[4], lwd = 1)
  if ("Fixed_tau_1"    %in% names(df)) lines(t, df$Fixed_tau_1,    col = cols[5], lwd = 1)
  if (identical(legend_corner, "bottomright")) {
    legend("bottomright", legend = method_labels, col = cols,
           lty = 1, lwd = c(2, 1.5, 1, 1, 1), bty = "n", cex = 0.85,
           inset = 0.02)
  } else {
    usr <- par("usr")
    legend(usr[1] + 0.03 * (usr[2] - usr[1]),
           usr[4] - 0.035 * (usr[4] - usr[3]),
           legend = method_labels, col = cols,
           lty = 1, lwd = c(2, 1.5, 1, 1, 1), bty = "n", cex = 0.85,
           xjust = 0, yjust = 1)
  }
  dev.off()
}

plot_mean_P2 <- function() {
  df <- read.csv(file.path(data_dir, "benchmark_P2_mean_trajectories_N1000.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  trajectory_plot(df, "Ensemble mean of P2",
                  file.path(out_dir, "fig_mean_P2.pdf"))
}

plot_sd_P2 <- function() {
  df <- read.csv(file.path(data_dir, "benchmark_P2_sd_trajectories_N1000.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  trajectory_plot(df, "Ensemble SD of P2",
                  file.path(out_dir, "fig_sd_P2.pdf"),
                  legend_corner = "bottomright")
}

# supervisor asked for a Pareto-style runtime vs accuracy plot in the
# meeting around mid-March. This is it.
plot_tradeoff <- function() {
  df <- summary_df
  lab <- ifelse(df$method == "SSA", "SSA",
         ifelse(startsWith(df$method, "Adaptive"), "Adaptive",
                paste0("Fixed tau=", df$tau)))
  rt    <- df$mean_runtime_pipeline_s
  rmm   <- df$RMSE_mean_P2
  rmsd  <- df$RMSE_sd_P2

  # two panels: runtime vs RMSE_mean on the left, runtime vs RMSE_sd on the right
  pdf(file.path(out_dir, "fig_tradeoff.pdf"), width = 9.5, height = 4)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  plot(rt, rmm, pch = 19, col = cols, cex = 1.5,
       xlim = c(0, max(rt) * 1.15),
       ylim = c(-1, max(rmm) * 1.15),
       xlab = "Mean runtime per replicate (s)",
       ylab = "RMSE of ensemble mean P2",
       main = "(a) RMSE on ensemble mean")
  text(rt, rmm, labels = lab, pos = 4, cex = 0.8)

  plot(rt, rmsd, pch = 19, col = cols, cex = 1.5,
       xlim = c(0, max(rt) * 1.15),
       ylim = c(-0.5, max(rmsd) * 1.15),
       xlab = "Mean runtime per replicate (s)",
       ylab = "RMSE of ensemble SD P2",
       main = "(b) RMSE on ensemble SD")
  text(rt, rmsd, labels = lab, pos = 4, cex = 0.8)
  dev.off()
}

cat("Reading from", data_dir, "\n")
cat("Writing to  ", out_dir, "\n")
plot_runtime();   cat("  fig_runtime.pdf\n")
plot_safeguard(); cat("  fig_safeguard.pdf\n")
plot_mean_P2();   cat("  fig_mean_P2.pdf\n")
plot_sd_P2();     cat("  fig_sd_P2.pdf\n")
plot_tradeoff();  cat("  fig_tradeoff.pdf\n")
cat("Done.\n")
