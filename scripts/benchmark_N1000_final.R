# N=1000 benchmark for chapter 5 (SSA, adaptive tau, fixed tau).
# Run from the benchmark/ folder:
#   Rscript scripts/benchmark_N1000_final.R               # full N=1000 -> outputs/
#   Rscript scripts/benchmark_N1000_final.R 2 outputs_test  # smoke test, 2 reps
# Writes per-method partial_*.rds first then merges to CSV at the end.

# force single thread BLAS so timing is fair
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

options(repos = c(CRAN = "https://cloud.r-project.org/"))

if (!requireNamespace("smfsb", quietly = TRUE)) install.packages("smfsb")
if (!requireNamespace("adaptivetau", quietly = TRUE)) install.packages("adaptivetau")
suppressPackageStartupMessages({
  library(smfsb)
  library(adaptivetau)
})

# use master_seed=42 same as project plan
master_seed <- 42L
N_rep       <- 1000L

# allow quick smoke test: first arg = N_rep, second arg = output dir.
# default still goes to outputs/ so the archived run is reproducible.
output_dir <- "outputs"
.cli <- commandArgs(trailingOnly = TRUE)
if (length(.cli) >= 1 && nzchar(.cli[1])) {
  v <- suppressWarnings(as.integer(.cli[1]))
  if (!is.na(v) && v >= 1) N_rep <- v
}
if (length(.cli) >= 2 && nzchar(.cli[2])) output_dir <- .cli[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("N_rep =", N_rep, "  output_dir =", output_dir, "\n")
t_grid      <- seq(0, 500, by = 1)
n_grid      <- length(t_grid)
tt          <- 500
dt          <- 1
adap_rrc    <- 0.15

x0 <- c(Gene = 10, P2Gene = 0, Rna = 0, P = 0, P2 = 0)
th <- c(k1 = 1, k1r = 10, k2 = 0.01, k3 = 10,
        k4 = 1, k4r = 1, k5 = 0.1, k6 = 0.01)

set.seed(master_seed)
seeds <- sample.int(.Machine$integer.max, N_rep)

init_vec   <- as.numeric(x0)
rrc_vec    <- rep(adap_rrc, length(x0))
p2_col_idx <- which(names(x0) == "P2") + 1L

# auto-regulatory model (smfsb form)
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

# safeguarded fixed-step tau leap (same as fixed_tau_safeguarded.R)
StepFixedTauLeap <- function(N, tau, safeguard_env,
                             max_shrink = 20L, min_tau = 1e-6) {
  S    <- t(N$Post - N$Pre)
  v    <- ncol(S)
  step_ssa <- StepGillespie(N)
  function(x0, t0, deltat, ...) {
    x <- x0
    t_remaining <- deltat
    while (t_remaining > 1e-10) {
      dt_sub <- min(tau, t_remaining)
      tau_try <- dt_sub
      accepted <- FALSE
      triggered <- FALSE
      for (k in seq_len(max_shrink)) {
        h <- pmax(N$h(x, t0 + (deltat - t_remaining), ...), 0)
        h[!is.finite(h)] <- 0
        K <- rpois(v, lambda = h * tau_try)
        x_new <- x + as.vector(S %*% K)
        if (!any(x_new < 0)) {
          x <- x_new
          t_remaining <- t_remaining - tau_try
          accepted <- TRUE
          break
        }
        # only count the first shrink in each sub-step
        if (!triggered) {
          safeguard_env$count <- safeguard_env$count + 1L
          triggered <- TRUE
        }
        tau_try <- tau_try / 2
        if (tau_try < min_tau) break
      }
      if (!accepted) {
        # fall back to one SSA step if shrink can not save it
        x <- step_ssa(x0 = x, t0 = t0 + (deltat - t_remaining),
                      deltat = dt_sub, ...)
        t_remaining <- t_remaining - dt_sub
      }
    }
    x
  }
}

# adaptivetau needs a different transition matrix layout
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

align_to_grid <- function(t_irreg, y_irreg, t_grid) {
  if (length(unique(t_irreg)) < 2) return(rep(NA_real_, length(t_grid)))
  approx(t_irreg, y_irreg, xout = t_grid,
         method = "constant", f = 0, rule = 2)$y
}

stepSSA <- StepGillespie(Auto)

run_ssa <- function(seed_i) {
  set.seed(seed_i)
  t0 <- proc.time()
  out <- simTs(x0 = Auto$M, t0 = 0, tt = tt, dt = dt,
               stepFun = stepSSA, th = th)
  P2 <- as.matrix(out)[, "P2"]
  t1 <- proc.time()
  list(P2 = P2,
       runtime_pipeline_s = unname((t1 - t0)["elapsed"]),
       runtime_core_s     = unname((t1 - t0)["elapsed"]),
       runtime_align_s    = 0,
       success = !any(is.na(P2)),
       safeguard_triggers = NA_integer_,
       n_events = NA_integer_)
}

run_adaptive <- function(seed_i) {
  set.seed(seed_i)
  t0 <- proc.time()
  raw <- ssa.adaptivetau(init.values = init_vec, transitions = nu,
                         rateFunc = rateFunc, params = th, tf = tt,
                         relratechange = rrc_vec)
  t_core <- proc.time()
  ok <- is.matrix(raw) && nrow(raw) >= 2 && raw[1, 1] == 0 &&
    ncol(raw) >= p2_col_idx
  if (ok) {
    P2 <- align_to_grid(raw[, 1], raw[, p2_col_idx], t_grid)
    n_ev <- nrow(raw)
  } else {
    P2 <- rep(NA_real_, n_grid)
    n_ev <- NA_integer_
  }
  t1 <- proc.time()
  list(P2 = P2,
       runtime_pipeline_s = unname((t1 - t0)["elapsed"]),
       runtime_core_s     = unname((t_core - t0)["elapsed"]),
       runtime_align_s    = unname((t1 - t_core)["elapsed"]),
       success = ok && !any(is.na(P2)),
       safeguard_triggers = NA_integer_,
       n_events = n_ev)
}

run_fixed <- function(seed_i, tau_j) {
  sg_env <- new.env(parent = emptyenv())
  sg_env$count <- 0L
  step_fn <- StepFixedTauLeap(Auto, tau = tau_j, safeguard_env = sg_env)
  set.seed(seed_i)
  t0 <- proc.time()
  out <- simTs(x0 = Auto$M, t0 = 0, tt = tt, dt = dt,
               stepFun = step_fn, th = th)
  P2 <- as.matrix(out)[, "P2"]
  t1 <- proc.time()
  list(P2 = P2,
       runtime_pipeline_s = unname((t1 - t0)["elapsed"]),
       runtime_core_s     = unname((t1 - t0)["elapsed"]),
       runtime_align_s    = 0,
       success = !any(is.na(P2)),
       safeguard_triggers = sg_env$count,
       n_events = NA_integer_)
}

method_order <- c("ssa", "adaptive", "fixed_0.25", "fixed_0.5", "fixed_1.0")
cat("N =", N_rep, " master_seed =", master_seed, "\n")
cat("outputs ->", output_dir, "\n\n")

# this loop is slow but easier to read than parallel
for (method_arg in method_order) {
  cat("--------\nmethod =", method_arg, "\n")
  P2_mat   <- matrix(NA_real_, nrow = n_grid, ncol = N_rep)
  rt_pipe  <- numeric(N_rep)
  rt_core  <- numeric(N_rep)
  rt_align <- numeric(N_rep)
  succ     <- logical(N_rep)
  sg_vec   <- integer(N_rep)
  n_ev_vec <- integer(N_rep)

  t_meth0 <- proc.time()[3]
  for (i in seq_len(N_rep)) {
    if (method_arg == "ssa") {
      r <- run_ssa(seeds[i])
    } else if (method_arg == "adaptive") {
      r <- run_adaptive(seeds[i])
    } else if (method_arg == "fixed_0.25") {
      r <- run_fixed(seeds[i], 0.25)
    } else if (method_arg == "fixed_0.5") {
      r <- run_fixed(seeds[i], 0.5)
    } else {
      r <- run_fixed(seeds[i], 1.0)
    }
    P2_mat[, i] <- r$P2
    rt_pipe[i]  <- r$runtime_pipeline_s
    rt_core[i]  <- r$runtime_core_s
    rt_align[i] <- r$runtime_align_s
    succ[i]     <- r$success
    sg_vec[i]   <- if (is.na(r$safeguard_triggers)) NA_integer_
                   else as.integer(r$safeguard_triggers)
    n_ev_vec[i] <- if (is.na(r$n_events)) NA_integer_
                   else as.integer(r$n_events)
    # print every 25 reps so we know it is not stuck
    if (i %% 25 == 0L || i == N_rep) {
      el <- proc.time()[3] - t_meth0
      cat(sprintf("  %4d/%d  pipe=%.3fs  total=%.0fs\n",
                  i, N_rep, rt_pipe[i], el))
      flush.console()
    }
  }

  tau_or_rrc <- switch(method_arg,
                       ssa = "-", adaptive = "uniform rrc=0.15",
                       fixed_0.25 = "0.25", fixed_0.5 = "0.5",
                       fixed_1.0 = "1")
  rep_df <- data.frame(
    method = method_arg,
    tau_or_rrc = tau_or_rrc,
    replicate_id = seq_len(N_rep),
    seed = seeds,
    runtime_pipeline_s = rt_pipe,
    runtime_core_s     = rt_core,
    runtime_align_s    = rt_align,
    success = succ,
    safeguard_triggers = sg_vec,
    n_events = n_ev_vec,
    stringsAsFactors = FALSE
  )

  partial <- list(
    method_key      = method_arg,
    master_seed     = master_seed,
    seeds           = seeds,
    t_grid          = t_grid,
    P2_matrix       = P2_mat,
    rep_df          = rep_df,
    R_version       = R.version.string,
    pkg_smfsb       = as.character(packageVersion("smfsb")),
    pkg_adaptivetau = as.character(packageVersion("adaptivetau")),
    sys_info        = as.list(Sys.info()),
    started_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  saveRDS(partial,
          file.path(output_dir, sprintf("partial_%s.rds", method_arg)))
  rm(partial, P2_mat, rep_df, rt_pipe, rt_core, rt_align,
     succ, sg_vec, n_ev_vec, t_meth0)
  gc(full = TRUE)
  cat("Saved partial_", method_arg, ".rds\n", sep = "")

  # 5-minute cool-down between methods so each method starts from a
  # comparable thermal/cache state. The cool-down is part of the
  # reproducibility protocol cited in Appendix A.6 and Section 4.4.1.
  if (!identical(method_arg, tail(method_order, 1L))) {
    cat(sprintf("  cool-down 5 min before next method (Sys.sleep(300))...\n"))
    flush.console()
    Sys.sleep(300)
    cat("  cool-down done, starting next method\n")
    flush.console()
  }
}

# --- merge step ---

required <- c("partial_ssa.rds", "partial_adaptive.rds",
              "partial_fixed_0.25.rds", "partial_fixed_0.5.rds",
              "partial_fixed_1.0.rds")
partials <- lapply(required, function(f) readRDS(file.path(output_dir, f)))
names(partials) <- gsub("^partial_|\\.rds$", "", required)

t_grid <- partials[[1]]$t_grid
N_rep  <- ncol(partials[[1]]$P2_matrix)

mean_ssa <- rowMeans(partials$ssa$P2_matrix, na.rm = TRUE)
sd_ssa   <- apply(partials$ssa$P2_matrix, 1, sd, na.rm = TRUE)

rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))

method_label <- function(k) {
  switch(k, ssa = "SSA", adaptive = "Adaptive tau", "Fixed tau")
}
tau_label <- function(k) {
  switch(k, ssa = "-", adaptive = "uniform rrc=0.15",
         fixed_0.25 = "0.25", fixed_0.5 = "0.5", fixed_1.0 = "1")
}

summary_rows <- list()
mean_traj <- data.frame(time = t_grid)
sd_traj   <- data.frame(time = t_grid)
all_rep   <- list()
all_sg    <- list()

for (k in c("ssa", "adaptive", "fixed_0.25", "fixed_0.5", "fixed_1.0")) {
  P <- partials[[k]]$P2_matrix
  m <- rowMeans(P, na.rm = TRUE)
  s <- apply(P, 1, sd, na.rm = TRUE)

  col_name <- switch(k,
                     ssa = "SSA", adaptive = "Adaptive",
                     fixed_0.25 = "Fixed_tau_0.25",
                     fixed_0.5 = "Fixed_tau_0.5",
                     fixed_1.0 = "Fixed_tau_1")
  mean_traj[[col_name]] <- m
  sd_traj[[col_name]]   <- s

  rdf <- partials[[k]]$rep_df

  summary_rows[[k]] <- data.frame(
    method = method_label(k),
    tau    = tau_label(k),
    n_rep  = N_rep,
    mean_runtime_pipeline_s   = mean(rdf$runtime_pipeline_s),
    median_runtime_pipeline_s = median(rdf$runtime_pipeline_s),
    sd_runtime_pipeline_s     = sd(rdf$runtime_pipeline_s),
    min_pipeline_s            = min(rdf$runtime_pipeline_s),
    max_pipeline_s            = max(rdf$runtime_pipeline_s),
    mean_runtime_core_s       = mean(rdf$runtime_core_s),
    mean_runtime_align_s      = mean(rdf$runtime_align_s),
    RMSE_mean_P2 = if (k == "ssa") 0 else rmse(m, mean_ssa),
    RMSE_sd_P2   = if (k == "ssa") 0 else rmse(s, sd_ssa),
    mean_safeguard_triggers = if (startsWith(k, "fixed"))
                                mean(rdf$safeguard_triggers, na.rm = TRUE)
                              else NA_real_,
    mean_n_events = if (k == "adaptive") mean(rdf$n_events, na.rm = TRUE)
                    else NA_real_,
    n_failures = sum(!rdf$success),
    stringsAsFactors = FALSE
  )
  all_rep[[k]] <- rdf
  if (startsWith(k, "fixed")) {
    all_sg[[k]] <- data.frame(
      tau = as.numeric(tau_label(k)),
      replicate_id = rdf$replicate_id,
      safeguard_triggers = rdf$safeguard_triggers
    )
  }
}

# tag the file name with N so a smoke run can not overwrite the N=1000 archive
n_tag <- sprintf("N%d", N_rep)
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df,
          file.path(output_dir, sprintf("benchmark_summary_%s.csv", n_tag)),
          row.names = FALSE)
write.csv(mean_traj,
          file.path(output_dir, sprintf("benchmark_P2_mean_trajectories_%s.csv", n_tag)),
          row.names = FALSE)
write.csv(sd_traj,
          file.path(output_dir, sprintf("benchmark_P2_sd_trajectories_%s.csv", n_tag)),
          row.names = FALSE)

rep_all <- do.call(rbind, all_rep)
write.csv(rep_all,
          file.path(output_dir, sprintf("replicate_level_runtime_%s.csv", n_tag)),
          row.names = FALSE)
sg_all <- do.call(rbind, all_sg)
write.csv(sg_all,
          file.path(output_dir, sprintf("replicate_level_safeguard_%s.csv", n_tag)),
          row.names = FALSE)

# metadata file. The appendix reproducibility table cites these rows by
# name, so keep BLAS threads, fairness_fixes, timing_convention and
# reference_for_RMSE in here even though we do not read them back.
meta <- data.frame(
  item = c("R_version", "pkg_smfsb", "pkg_adaptivetau",
           "sys_info_sysname", "sys_info_release",
           "sys_info_nodename", "sys_info_machine",
           "detectCores", "blas_threads",
           "master_seed", "n_rep", "adaptive_relratechange",
           "reference_for_RMSE", "timing_convention",
           "fairness_fixes", "merged_at"),
  value = c(partials[[1]]$R_version,
            partials[[1]]$pkg_smfsb,
            partials[[1]]$pkg_adaptivetau,
            as.character(partials[[1]]$sys_info$sysname),
            as.character(partials[[1]]$sys_info$release),
            as.character(partials[[1]]$sys_info$nodename),
            as.character(partials[[1]]$sys_info$machine),
            as.character(parallel::detectCores()),
            "OMP=1,MKL=1,OPENBLAS=1,VECLIB=1",
            as.character(partials[[1]]$master_seed),
            as.character(N_rep),
            "rep(0.15, 5) all species = uniform 0.15",
            "SSA ensemble mean / sd from the same per-replicate seeds",
            "pipeline (x0 -> P2 vector on common grid, includes align_to_grid)",
            "single-process cold-start sequential; 5-min cool-down between methods",
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
meta_file <- if (N_rep == 1000L) "run_metadata_N1000_v3.csv" else sprintf("run_metadata_%s.csv", n_tag)
write.csv(meta, file.path(output_dir, meta_file), row.names = FALSE)

# keep partial_*.rds on disk; chapter5_postprocess.R needs the per-method
# P2 matrices for the across-replicate paired bootstrap.
cat("Merged into benchmark_*.csv. Partials kept for chapter5_postprocess.R.\n")
print(summary_df, row.names = FALSE, digits = 4)
