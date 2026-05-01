# One representative SSA trajectory of the auto-regulatory model for the
# chapter 3 figure. Writes figures/fig_sample_trajectory.pdf. Uses the
# first seed from the master_seed=42 vector so the displayed path is the
# same one across reruns.

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

suppressPackageStartupMessages(library(smfsb))

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

x0 <- c(Gene = 10, P2Gene = 0, Rna = 0, P = 0, P2 = 0)
th <- c(k1 = 1, k1r = 10, k2 = 0.01, k3 = 10,
        k4 = 1, k4r = 1, k5 = 0.1, k6 = 0.01)

Auto <- list()
Auto$M <- x0
Auto$Pre <- matrix(
  c(1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0),
  ncol = 5, byrow = TRUE)
colnames(Auto$Pre) <- names(Auto$M)
Auto$Post <- matrix(
  c(0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  ncol = 5, byrow = TRUE)
colnames(Auto$Post) <- names(Auto$M)
local({
  th_local <- th
  Auto$h <<- function(x, t, th = th_local) {
    x <- setNames(as.numeric(x), names(Auto$M))
    with(as.list(c(x, th)),
         c(k1 * Gene * P2, k1r * P2Gene, k2 * Gene, k3 * Rna,
           k4 * P * (P - 1) / 2, k4r * P2, k5 * Rna, k6 * P))
  }
})

set.seed(42L)
seeds <- sample.int(.Machine$integer.max, 5L)
set.seed(seeds[1])

stepSSA <- StepGillespie(Auto)
out <- simTs(x0 = Auto$M, t0 = 0, tt = 500, dt = 0.5,
             stepFun = stepSSA, th = th)
M <- as.matrix(out)
t_grid <- as.numeric(time(out))

pdf(file.path(out_dir, "fig_sample_trajectory.pdf"),
    width = 6.5, height = 4.2)
par(mar = c(4, 4, 1, 1), las = 1)
plot(t_grid, M[, "P2"], type = "l", col = "darkgreen", lwd = 1.4,
     xlab = "Time (s)", ylab = "Copy number",
     ylim = c(0, max(M[, c("P2", "P", "Rna")]) * 1.05))
lines(t_grid, M[, "P"],   col = "blue",   lwd = 1, lty = 2)
lines(t_grid, M[, "Rna"], col = "orange", lwd = 1, lty = 3)
lines(t_grid, M[, "P2Gene"] * 20, col = "red", lwd = 1, lty = 1)
usr <- par("usr")
legend(usr[2] - 0.03 * (usr[2] - usr[1]),
       usr[3] + 0.50 * (usr[4] - usr[3]),
       legend = c(expression(P[2]), "P (monomer)", "R (mRNA)",
                  expression(G[b] %*% 20 ~~ "(promoter bound)")),
       col = c("darkgreen", "blue", "orange", "red"),
       lty = c(1, 2, 3, 1), lwd = c(1.4, 1, 1, 1), bty = "n",
       cex = 0.85, xjust = 1, yjust = 0.5)
dev.off()
cat("Wrote", file.path(out_dir, "fig_sample_trajectory.pdf"), "\n")
