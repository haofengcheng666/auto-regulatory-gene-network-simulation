# Terminal-time ECDF and KS distance for the chapter 5
# distribution-level check. Reads outputs/partial_*.rds, takes the
# P_2(t=500) column from each per-method matrix, and writes
#   figures/fig_terminal_ecdf.pdf
#   outputs/ks_terminal_N1000.csv  (one row per method)

data_dir <- "outputs"
out_dir  <- "figures"
.cli <- commandArgs(trailingOnly = TRUE)
if (length(.cli) >= 1 && nzchar(.cli[1])) data_dir <- .cli[1]
if (length(.cli) >= 2 && nzchar(.cli[2])) out_dir  <- .cli[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required <- c("partial_ssa.rds", "partial_adaptive.rds",
              "partial_fixed_0.25.rds", "partial_fixed_0.5.rds",
              "partial_fixed_1.0.rds")
missing <- required[!file.exists(file.path(data_dir, required))]
if (length(missing) > 0L) {
  stop("Missing partials in ", data_dir, ": ",
       paste(missing, collapse = ", "),
       ". Run benchmark_N1000_final.R first.")
}

partials <- lapply(required, function(f) readRDS(file.path(data_dir, f)))
names(partials) <- gsub("^partial_|\\.rds$", "", required)

t_grid <- partials[[1]]$t_grid
i500   <- which(t_grid == 500)
stopifnot(length(i500) == 1L)

# pull terminal-time P_2 values, one length-N vector per method
P2_500 <- lapply(partials, function(p) as.numeric(p$P2_matrix[i500, ]))

method_keys   <- c("ssa", "adaptive", "fixed_0.25", "fixed_0.5", "fixed_1.0")
method_labels <- c("SSA", "Adaptive tau",
                   "Fixed tau=0.25", "Fixed tau=0.5", "Fixed tau=1.0")
cols <- c("black", "blue", "orange", "red", "darkgreen")

# ECDF panel
pdf(file.path(out_dir, "fig_terminal_ecdf.pdf"), width = 6.5, height = 4)
par(mar = c(4.2, 4.2, 1, 1), las = 1)
xrange <- range(unlist(P2_500))
plot(NA, xlim = xrange, ylim = c(0, 1), xlab = expression(P[2](t == 500)),
     ylab = "Empirical CDF")
for (i in seq_along(method_keys)) {
  lines(ecdf(P2_500[[method_keys[i]]]), col = cols[i],
        lwd = if (method_keys[i] == "ssa") 2 else 1.4,
        do.points = FALSE, verticals = TRUE)
}
legend("bottomright", legend = method_labels, col = cols,
       lty = 1, lwd = c(2, 1.4, 1.4, 1.4, 1.4), bty = "n", cex = 0.85)
dev.off()
cat("Wrote", file.path(out_dir, "fig_terminal_ecdf.pdf"), "\n")

# KS distances vs SSA reference
ssa_vec <- P2_500[["ssa"]]
ks_rows <- list()
ks_rows[["ssa"]] <- data.frame(
  method = "SSA", tau_or_rrc = "-", n_rep = length(ssa_vec),
  mean_P2_500 = mean(ssa_vec), sd_P2_500 = sd(ssa_vec),
  KS_vs_SSA = 0, KS_pvalue = NA_real_,
  stringsAsFactors = FALSE
)
for (k in c("adaptive", "fixed_0.25", "fixed_0.5", "fixed_1.0")) {
  v <- P2_500[[k]]
  kk <- suppressWarnings(ks.test(v, ssa_vec))
  ks_rows[[k]] <- data.frame(
    method = switch(k, adaptive = "Adaptive tau",
                    fixed_0.25 = "Fixed tau=0.25",
                    fixed_0.5  = "Fixed tau=0.5",
                    fixed_1.0  = "Fixed tau=1.0"),
    tau_or_rrc = switch(k, adaptive = "0.15", fixed_0.25 = "0.25",
                        fixed_0.5 = "0.5", fixed_1.0 = "1"),
    n_rep = length(v),
    mean_P2_500 = mean(v), sd_P2_500 = sd(v),
    KS_vs_SSA = unname(kk$statistic),
    KS_pvalue = unname(kk$p.value),
    stringsAsFactors = FALSE
  )
}
out_csv <- file.path(data_dir, sprintf("ks_terminal_N%d.csv", length(ssa_vec)))
ks_df <- do.call(rbind, ks_rows)
write.csv(ks_df, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")

cat("\n----- terminal-time summary -----\n")
print(ks_df, row.names = FALSE, digits = 4)
