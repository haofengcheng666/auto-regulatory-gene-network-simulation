# Parallel adaptive-tau tolerance pilot for appendix A.5.
#
# Drop-in faster version of pilot_sweep.R. Runs adaptive tau-leaping at
# rrc in {0.05, 0.10, 0.15, 0.20} and an SSA reference, all on the same
# auto-regulatory model and the same first N_pilot seeds as the main
# N=1000 benchmark, but distributes replicates across n_workers PSOCK
# workers.
#
# Each setting's full per-replicate output is saved to a single RDS file
# as soon as it completes, so a crash in any later step never costs more
# than the most recent setting.
#
# IMPORTANT: pilot runtime in parallel mode is the per-replicate R-level
# core-time measured inside each worker (proc.time() on a single replicate
# inside an isolated R process). It is NOT the wall-clock pipeline time
# that appears in Section 5.1 of the main benchmark, which used a strict
# single-process cold-start protocol. The pilot's primary purpose is the
# RMSE comparison across rrc; the runtime column is reported as a
# coarse-grained sanity check only.
#
# Usage:
#   Rscript scripts/pilot_sweep_parallel.R                  # N=100, 8 workers
#   Rscript scripts/pilot_sweep_parallel.R --N=1000         # full size
#   Rscript scripts/pilot_sweep_parallel.R --N=20 --workers=4 outputs_smoke

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

suppressPackageStartupMessages({
  library(smfsb)
  library(adaptivetau)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
N_pilot <- 100L
n_workers <- min(8L, parallel::detectCores(logical = FALSE))
output_dir <- "outputs"
for (a in args) {
  if (grepl("^--N=", a)) {
    N_pilot <- as.integer(sub("^--N=", "", a))
  } else if (grepl("^--workers=", a)) {
    n_workers <- as.integer(sub("^--workers=", "", a))
  } else if (nzchar(a) && !startsWith(a, "--")) {
    output_dir <- a
  }
}
cat(sprintf("pilot_sweep_parallel.R  N_pilot = %d  n_workers = %d  output_dir = %s\n",
            N_pilot, n_workers, output_dir))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

master_seed <- 42L
t_grid <- seq(0, 500, by = 1)
n_grid <- length(t_grid)
tt <- 500
dt <- 1

x0 <- c(Gene = 10, P2Gene = 0, Rna = 0, P = 0, P2 = 0)
th <- c(k1 = 1, k1r = 10, k2 = 0.01, k3 = 10,
        k4 = 1, k4r = 1, k5 = 0.1, k6 = 0.01)

set.seed(master_seed)
seeds <- sample.int(.Machine$integer.max, N_pilot)

init_vec <- as.numeric(x0)
p2_col_idx <- which(names(x0) == "P2") + 1L

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
stepSSA <- StepGillespie(Auto)

nu <- cbind(c(-1, 1, 0, 0, -1), c(1, -1, 0, 0, 1), c(0, 0, 1, 0, 0),
            c(0, 0, 0, 1, 0), c(0, 0, 0, -2, 1), c(0, 0, 0, 2, -1),
            c(0, 0, -1, 0, 0), c(0, 0, 0, -1, 0))
rownames(nu) <- names(x0)

rateFunc <- function(x, params, t) {
  xx <- pmax(setNames(as.numeric(x), names(x0)), 0)
  with(as.list(c(xx, params)),
       c(k1 * Gene * P2, k1r * P2Gene, k2 * Gene, k3 * Rna,
         k4 * P * (P - 1) / 2, k4r * P2, k5 * Rna, k6 * P))
}

run_ssa <- function(seed_i) {
  set.seed(seed_i)
  t0 <- proc.time()
  out <- simTs(x0 = Auto$M, t0 = 0, tt = tt, dt = dt,
               stepFun = stepSSA, th = th)
  P2 <- as.matrix(out)[, "P2"]
  list(P2 = P2,
       runtime_pipeline_s = unname((proc.time() - t0)["elapsed"]),
       n_events = NA_integer_)
}

run_adaptive <- function(seed_i, rrc) {
  rrc_vec <- rep(rrc, length(x0))
  set.seed(seed_i)
  t0 <- proc.time()
  raw <- ssa.adaptivetau(init.values = init_vec, transitions = nu,
                         rateFunc = rateFunc, params = th, tf = tt,
                         relratechange = rrc_vec)
  ok <- is.matrix(raw) && nrow(raw) >= 2 && raw[1, 1] == 0 &&
    ncol(raw) >= p2_col_idx
  if (ok) {
    P2 <- approx(raw[, 1], raw[, p2_col_idx], xout = t_grid,
                 method = "constant", f = 0, rule = 2)$y
    n_ev <- nrow(raw)
  } else {
    P2 <- rep(NA_real_, n_grid)
    n_ev <- NA_integer_
  }
  list(P2 = P2,
       runtime_pipeline_s = unname((proc.time() - t0)["elapsed"]),
       n_events = n_ev)
}

# ---- cluster setup ---------------------------------------------------------
cat(sprintf("Spinning up %d PSOCK workers ...\n", n_workers))
cl <- makeCluster(n_workers, type = "PSOCK")
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  suppressPackageStartupMessages({
    library(smfsb)
    library(adaptivetau)
  })
  TRUE
}))
clusterExport(cl, varlist = c(
  "Auto", "stepSSA", "nu", "rateFunc",
  "init_vec", "p2_col_idx", "x0", "th",
  "tt", "dt", "t_grid", "n_grid",
  "run_ssa", "run_adaptive"
), envir = environment())
cat("Cluster ready.\n")

# run_block_parallel passes the run-function and any extra args (e.g. rrc)
# directly through parLapplyLB's `...`, which is the cleanest way to avoid
# closure-serialisation pitfalls across PSOCK workers.
run_block_parallel <- function(label, fn, ...) {
  t0 <- proc.time()[3]
  results <- parLapplyLB(cl, seeds, fn, ...)
  P2_mat <- do.call(cbind, lapply(results, `[[`, "P2"))
  rt <- vapply(results, `[[`, numeric(1), "runtime_pipeline_s")
  ne <- vapply(results, function(r) {
    v <- r$n_events
    if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v)
  }, integer(1))
  wall <- proc.time()[3] - t0
  cat(sprintf("  %s  done  N=%d  wall=%.0fs  mean(core)=%.3fs  workers=%d  speedup~%.2fx\n",
              label, N_pilot, wall, mean(rt), n_workers,
              (mean(rt) * N_pilot) / wall))
  flush.console()
  list(P2 = P2_mat, runtime_s = rt, n_events = ne)
}

save_rds_safe <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
  cat("Saved", path, "\n")
}

# ---- SSA reference ---------------------------------------------------------
ssa_path <- file.path(output_dir, sprintf("pilot_partial_ssa_N%d.rds", N_pilot))
if (file.exists(ssa_path)) {
  cat("Skipping SSA: already at", ssa_path, "\n")
} else {
  cat("---- SSA reference (N =", N_pilot, ") ----\n")
  ssa_block <- run_block_parallel("ssa", run_ssa)
  save_rds_safe(list(method = "ssa", N_pilot = N_pilot, n_workers = n_workers,
                     master_seed = master_seed, seeds = seeds,
                     t_grid = t_grid,
                     P2 = ssa_block$P2,
                     runtime_s = ssa_block$runtime_s,
                     n_events = ssa_block$n_events),
                ssa_path)
}

# ---- adaptive at four rrc settings ----------------------------------------
rrc_grid <- c(0.05, 0.10, 0.15, 0.20)
for (rrc in rrc_grid) {
  rrc_path <- file.path(output_dir,
                        sprintf("pilot_partial_rrc_%.2f_N%d.rds",
                                rrc, N_pilot))
  if (file.exists(rrc_path)) {
    cat("Skipping rrc =", rrc, ": already at", rrc_path, "\n")
    next
  }
  cat(sprintf("---- adaptive  rrc = %.2f  (N = %d) ----\n", rrc, N_pilot))
  block <- run_block_parallel(sprintf("rrc=%.2f", rrc),
                              run_adaptive, rrc = rrc)
  save_rds_safe(list(method = "adaptive", rrc = rrc, N_pilot = N_pilot,
                     n_workers = n_workers,
                     master_seed = master_seed, seeds = seeds,
                     t_grid = t_grid,
                     P2 = block$P2,
                     runtime_s = block$runtime_s,
                     n_events = block$n_events),
                rrc_path)
}

cat("All simulation blocks saved. Run scripts/pilot_postprocess.R --N=",
    N_pilot, " to assemble the table.\n", sep = "")
