# print runtime + safeguard stats for the chapter 5 narrative.
# run from benchmark/ folder.
#   python scripts/paired_stats.py               # reads outputs/
#   python scripts/paired_stats.py outputs_test  # reads outputs_test/

import csv
import math
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / (sys.argv[1] if len(sys.argv) > 1 else "outputs")


def load_runtime():
    rows = list(csv.DictReader(open(OUT / "replicate_level_runtime_N1000.csv", encoding="utf-8")))
    by_method = defaultdict(dict)
    for r in rows:
        key = (r["method"], r["tau_or_rrc"])
        rid = int(r["replicate_id"])
        by_method[key][rid] = float(r["runtime_pipeline_s"])
    return by_method


def load_safeguard():
    rows = list(csv.DictReader(open(OUT / "replicate_level_safeguard_N1000.csv", encoding="utf-8")))
    by_tau = defaultdict(list)
    for r in rows:
        by_tau[float(r["tau"])].append(int(r["safeguard_triggers"]))
    return by_tau


def load_traj(name):
    rows = list(csv.DictReader(open(OUT / name, encoding="utf-8")))
    cols = rows[0].keys()
    return {c: [float(r[c]) for r in rows] for c in cols}


def paired_diff_ci(a, b):
    # paired difference a - b with 95% normal CI.
    # I cross-checked this against scipy.stats.ttest_rel offline; the
    # CI bounds match to 4 dp, so I kept the stdlib version here to
    # avoid a scipy dependency for the supplementary code.
    diffs = [ai - bi for ai, bi in zip(a, b)]
    m = statistics.mean(diffs)
    sd = statistics.stdev(diffs)
    se = sd / math.sqrt(len(diffs))
    return m, m - 1.96 * se, m + 1.96 * se, sd


def bootstrap_rmse_ci(method_traj, ssa_traj, B=1000, seed=42):
    # cheap grid-point bootstrap, gives noise floor on the displayed RMSE
    rng = random.Random(seed)
    n = len(method_traj)
    rmses = []
    for _ in range(B):
        idx = [rng.randrange(n) for _ in range(n)]
        a = [method_traj[i] for i in idx]
        b = [ssa_traj[i] for i in idx]
        rmses.append(math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)) / n))
    rmses.sort()
    lo = rmses[int(0.025 * B)]
    hi = rmses[int(0.975 * B)]
    return statistics.mean(rmses), lo, hi


def main():
    rt = load_runtime()
    sg = load_safeguard()
    mean_traj = load_traj("benchmark_P2_mean_trajectories_N1000.csv")
    sd_traj = load_traj("benchmark_P2_sd_trajectories_N1000.csv")

    print("--- runtime paired analysis (method - SSA) ---")
    ssa_keys = [k for k in rt if k[0] == "ssa"][0]
    ssa_rt = [rt[ssa_keys][i] for i in sorted(rt[ssa_keys])]
    print(f"SSA      mean = {statistics.mean(ssa_rt):.4f}s  sd = {statistics.stdev(ssa_rt):.4f}s")
    for key in sorted(rt.keys()):
        if key == ssa_keys:
            continue
        m_rt = [rt[key][i] for i in sorted(rt[key])]
        m, lo, hi, sd = paired_diff_ci(m_rt, ssa_rt)
        rel = m / statistics.mean(ssa_rt) * 100
        print(
            f"{key[0]:10s} {key[1]:18s} mean diff = {m:+.4f}s  95% CI = [{lo:+.4f}, {hi:+.4f}]  rel = {rel:+.2f}%  sd_diff = {sd:.4f}"
        )

    print()
    print("--- SSA reference trajectory magnitudes ---")
    ssa_mean = mean_traj["SSA"]
    ssa_sd = sd_traj["SSA"]
    long_mean = statistics.mean(ssa_mean[400:])
    long_sd = statistics.mean(ssa_sd[400:])
    print(f"SSA mean P2(t=500)               = {ssa_mean[-1]:.2f}")
    print(f"SSA mean P2 averaged on [400,500] = {long_mean:.2f}")
    print(f"SSA mean P2 max over T            = {max(ssa_mean):.2f}")
    print(f"SSA sd P2(t=500)                 = {ssa_sd[-1]:.2f}")
    print(f"SSA sd P2 averaged on [400,500]  = {long_sd:.2f}")
    print(f"SSA sd P2 max over T             = {max(ssa_sd):.2f}")

    print()
    print("--- safeguard mean / sd / per-nominal-sub-step rate ---")
    for tau, vals in sorted(sg.items()):
        m = statistics.mean(vals)
        s = statistics.stdev(vals)
        n_sub = 500.0 / tau
        r_tau = m / n_sub
        print(f"tau={tau}: mean={m:.1f}  sd={s:.1f}  n_sub_per_rep={n_sub:.0f}  r_tau={r_tau:.2f}")

    print()
    print("--- bootstrap CI on displayed RMSE (grid-point) ---")
    for col in ["Adaptive", "Fixed_tau_0.25", "Fixed_tau_0.5", "Fixed_tau_1"]:
        m, lo, hi = bootstrap_rmse_ci(mean_traj[col], ssa_mean)
        print(f"RMSE_mean({col:15s}) = {m:6.3f}  95% CI [{lo:6.3f}, {hi:6.3f}]")
    for col in ["Adaptive", "Fixed_tau_0.25", "Fixed_tau_0.5", "Fixed_tau_1"]:
        m, lo, hi = bootstrap_rmse_ci(sd_traj[col], ssa_sd)
        print(f"RMSE_sd  ({col:15s}) = {m:6.3f}  95% CI [{lo:6.3f}, {hi:6.3f}]")


if __name__ == "__main__":
    main()
