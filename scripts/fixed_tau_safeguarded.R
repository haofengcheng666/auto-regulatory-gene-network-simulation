# Reference implementation of safeguarded fixed-step tau leap.
# Same step logic as the N=1000 benchmark, just one replicate for plotting.

set.seed(1)

library(smfsb)

Auto <- list()
Auto$M <- c(Gene = 10, P2Gene = 0, Rna = 0, P = 0, P2 = 0)

Auto$Pre <- matrix(
  c(
    1, 0, 0, 0, 1,
    0, 1, 0, 0, 0,
    1, 0, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 2, 0,
    0, 0, 0, 0, 1,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0
  ),
  ncol = 5, byrow = TRUE
)
colnames(Auto$Pre) <- names(Auto$M)

Auto$Post <- matrix(
  c(
    0, 1, 0, 0, 0,
    1, 0, 0, 0, 1,
    1, 0, 1, 0, 0,
    0, 0, 1, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 2, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  ),
  ncol = 5, byrow = TRUE
)
colnames(Auto$Post) <- names(Auto$M)

th <- c(k1 = 1, k1r = 10, k2 = 0.01, k3 = 10,
        k4 = 1, k4r = 1, k5 = 0.1, k6 = 0.01)

Auto$h <- function(x, t, th = th) {
  with(as.list(c(x, th)),
       c(k1 * Gene * P2, k1r * P2Gene, k2 * Gene, k3 * Rna,
         k4 * P * (P - 1) / 2, k4r * P2, k5 * Rna, k6 * P))
}

# only count the first shrink in each sub-step
StepFixedTauLeap <- function(N, tau, safeguard_env,
                             max_shrink = 20L, min_tau = 1e-6) {
  S    <- t(N$Post - N$Pre)
  v    <- ncol(S)
  step_ssa <- StepGillespie(N)

  function(x0, t0, deltat, ...) {
    x           <- x0
    t_remaining <- deltat

    while (t_remaining > 1e-10) {
      dt_sub  <- min(tau, t_remaining)
      tau_try <- dt_sub
      accepted <- FALSE
      triggered <- FALSE

      for (k in seq_len(max_shrink)) {
        h <- N$h(x, t0 + (deltat - t_remaining), ...)
        h <- pmax(h, 0)
        h[!is.finite(h)] <- 0
        K <- rpois(v, lambda = h * tau_try)
        x_new <- x + as.vector(S %*% K)

        if (!any(x_new < 0)) {
          x          <- x_new
          t_remaining <- t_remaining - tau_try
          accepted   <- TRUE
          break
        }

        if (!triggered) {
          safeguard_env$count <- safeguard_env$count + 1L
          triggered           <- TRUE
        }
        tau_try <- tau_try / 2
        if (tau_try < min_tau) break
      }

      if (!accepted) {
        # SSA fallback if shrink can not save it
        x <- step_ssa(x0 = x, t0 = t0 + (deltat - t_remaining),
                      deltat = dt_sub, ...)
        t_remaining <- t_remaining - dt_sub
      }
    }

    x
  }
}

# quick visual check, not for the dissertation figures.
# one replicate at tau = 1 just so I can eyeball the species traces.
tau <- 1
sg_env <- new.env(parent = emptyenv())
sg_env$count <- 0L

stepTL <- StepFixedTauLeap(Auto, tau = tau, safeguard_env = sg_env)

tt <- 500
dt <- tau
outTL <- simTs(x0 = Auto$M, t0 = 0, tt = tt, dt = dt,
               stepFun = stepTL, th = th)

cat("Safeguard triggers (total):", sg_env$count, "\n")
plot(outTL)

stepSSA <- StepGillespie(Auto)
outSSA  <- simTs(x0 = Auto$M, t0 = 0, tt = tt, dt = dt,
                 stepFun = stepSSA, th = th)

matTL  <- as.matrix(outTL)
matSSA <- as.matrix(outSSA)
time <- time(outTL)

# compare three species, SSA solid vs tau leap dashed
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(time, matSSA[, "Rna"], type = "l", xlab = "Time", ylab = "Rna",
     main = "SSA (solid) vs fixed tau (dashed)")
lines(time, matTL[, "Rna"], lty = 2)
plot(time, matSSA[, "P"], type = "l", xlab = "Time", ylab = "P",
     main = "SSA (solid) vs fixed tau (dashed)")
lines(time, matTL[, "P"], lty = 2)
plot(time, matSSA[, "P2"], type = "l", xlab = "Time", ylab = "P2",
     main = "SSA (solid) vs fixed tau (dashed)")
lines(time, matTL[, "P2"], lty = 2)
par(mfrow = c(1, 1))
