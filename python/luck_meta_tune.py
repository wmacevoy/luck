"""Meta-tune the online luck tracker against a pool dataset.

Three-stage harness:

  1. Global fit on the pool (luck_fit.LuckFit) -> theta_seed, baseline rms_L.
     This is the "factory calibration" the online tracker starts from.

  2. Bootstrap stationary noise floors.  For each candidate window W, draw
     many independent length-W blocks from the model at theta_seed, run
     the EWMA-on-L update, and measure sd(ewma_mean), sd(ewma_d2) under
     a *correctly-fit, stationary* world.  These replace the closed-form
     1/sqrt(24 W) and 1/sqrt(360 W) used by LuckTrackerOnline so the
     trigger thresholds reflect the actual model class and finite M, not
     the asymptotic uniform.

  3. Synthetic-drift sweep.  Run the tracker on a known scale-step
     scenario for each (W, k_sigma) pair and score by composite

         score = lag_rms_in_theta + lambda_cost * refits_per_sample.

  Output: print the table; write the chosen settings to luck_config.json
  for LuckTrackerOnline to consume.

Run the demo block at the bottom for a self-contained pool of multivariate
Gaussian samples in df=12; swap that out to point the harness at any
saved data.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from luck_fit import LuckFit, mvt_family
from luck_tracker_online import LuckTrackerOnline


def global_fit(X_pool, theta0, n_sim=4000):
    df = X_pool.shape[1]
    mu_hat = X_pool.mean(axis=0)
    Sigma_hat = np.cov(X_pool.T, bias=False)
    log_pdf, sampler = mvt_family(mu_hat, Sigma_hat)
    lf = LuckFit(log_pdf, sampler, theta0=theta0, n_sim=n_sim)
    lf.add_samples(X_pool)
    lf.fit()
    return {
        "mu_hat": mu_hat, "Sigma_hat": Sigma_hat,
        "log_pdf": log_pdf, "sampler": sampler,
        "theta_seed": lf.theta.copy(),
        "rms_L_pool": lf.rms_error(lf.theta),
    }


def bootstrap_noise_floors(log_pdf, sampler, theta, W, M, n_trials=80, seed=1):
    """Empirical sd of EWMA(mean L) and EWMA((L - 1/2)^2) under stationarity."""
    rng = np.random.default_rng(seed)
    a = 1.0 / W
    ewma_means = np.empty(n_trials)
    ewma_d2s = np.empty(n_trials)
    X_sim = sampler(theta, M, rng)
    nlp_sim_sorted = np.sort(-log_pdf(theta, X_sim))
    for t in range(n_trials):
        X = sampler(theta, W, rng)
        nlp = -log_pdf(theta, X)
        idx_lo = np.searchsorted(nlp_sim_sorted, nlp, side="left")
        idx_hi = np.searchsorted(nlp_sim_sorted, nlp, side="right")
        L = 0.5 * (idx_lo + idx_hi) / M
        d = L - 0.5
        em = 0.5
        ed = 1.0 / 12.0
        for ell, dd in zip(L, d):
            em = (1 - a) * em + a * ell
            ed = (1 - a) * ed + a * dd * dd
        ewma_means[t] = em
        ewma_d2s[t] = ed
    return float(ewma_means.std(ddof=1)), float(ewma_d2s.std(ddof=1))


def make_drift_scenario(rng, mu_true, L_true, n_total, change_at, s_before, s_after):
    df = mu_true.size
    s = np.where(np.arange(n_total) < change_at, s_before, s_after)
    Z = rng.standard_normal((n_total, df))
    return mu_true + (s[:, None] * (Z @ L_true.T)), s


def run_tracker(log_pdf, sampler, theta_seed, df, X_drift, s_truth,
                W, k_sigma, sd_mean, sd_d2, sim_mult=4):
    tr = LuckTrackerOnline(log_pdf, sampler, theta0=theta_seed,
                           df=df, k_sigma=k_sigma, sim_mult=sim_mult)
    tr.W = W
    tr.alpha = 1.0 / W
    tr.refit_floor = W
    tr.burn_in = max(1, W // 20)
    tr.M = sim_mult * W
    tr._refresh_sim_buffer()

    def trigger_check(self_):
        if self_.t - self_.last_refit_t <= self_.burn_in:
            return None
        if abs(self_.ewma_mean - 0.5) > k_sigma * sd_mean:
            return "mean"
        if abs(self_.ewma_d2 - 1.0 / 12.0) > k_sigma * sd_d2:
            return "var"
        return None

    s_path = np.empty(len(X_drift))
    t0 = time.perf_counter()
    for i, x in enumerate(X_drift):
        tr.t += 1
        tr.window.append(x.copy())
        ell = tr._ell(x)
        d = ell - 0.5
        a_ = tr.alpha
        tr.ewma_mean = (1 - a_) * tr.ewma_mean + a_ * ell
        tr.ewma_d2 = (1 - a_) * tr.ewma_d2 + a_ * d * d
        tr.steps_since_refit += 1
        reason = trigger_check(tr)
        if reason or tr.steps_since_refit >= tr.refit_floor:
            tr._local_refit(reason or "floor")
        s_path[i] = tr.theta[0]
    elapsed = time.perf_counter() - t0

    settle_idx = int(0.6 * len(X_drift))
    s_err = s_path[settle_idx:] - s_truth[settle_idx:]
    return {
        "rms_s_error": float(np.sqrt((s_err ** 2).mean())),
        "refits": tr.refits,
        "elapsed": elapsed,
        "s_path_final": float(s_path[-1]),
    }


def main():
    rng = np.random.default_rng(20260430)
    df = 12
    A = rng.standard_normal((df, df))
    Sigma_true = A @ A.T / df + 0.5 * np.eye(df)
    mu_true = rng.standard_normal(df)
    L_true = np.linalg.cholesky(Sigma_true)

    n_pool = 4000
    X_pool = mu_true + (rng.standard_normal((n_pool, df)) @ L_true.T)

    print("=== stage 1: global fit on pool ===")
    g = global_fit(X_pool, theta0=[1.0, np.log(60)])
    print(f"  pool n={n_pool}, df={df}")
    print(f"  theta_seed = (s={g['theta_seed'][0]:.4f}, "
          f"nu={np.exp(g['theta_seed'][1]):.1f})")
    print(f"  rms_L on pool = {g['rms_L_pool']:.5f}\n")

    W_grid = [12 * df ** 2, 24 * df ** 2, 48 * df ** 2]
    k_grid = [2.0, 3.0, 4.0]
    sim_mult = 4

    print("=== stage 2: bootstrap stationary noise floors ===")
    print(f"  {'W':>6}   {'sd(mean) emp':>13}  {'sd(d^2) emp':>13}  "
          f"{'closed-form mean':>18}  {'closed-form d^2':>17}")
    floors = {}
    for W in W_grid:
        sdm, sdd = bootstrap_noise_floors(g["log_pdf"], g["sampler"],
                                          g["theta_seed"], W, M=sim_mult * W)
        floors[W] = (sdm, sdd)
        cf_m = np.sqrt((1 / 12) / (2 * W))
        cf_d = np.sqrt((1 / 180) / (2 * W))
        print(f"  {W:>6}   {sdm:>13.5f}  {sdd:>13.5f}  "
              f"{cf_m:>18.5f}  {cf_d:>17.5f}")

    print("\n=== stage 3: synthetic drift scoring (scale step 1.0 -> 1.4) ===")
    n_run = 6 * 24 * df ** 2
    X_drift, s_truth = make_drift_scenario(
        rng, mu_true, L_true, n_run, change_at=n_run // 3,
        s_before=1.0, s_after=1.4)

    lambda_cost = 0.005
    rows = []
    print(f"  scenario length n={n_run}, change at t={n_run // 3}")
    print(f"  {'W':>6}  {'k_sigma':>7}  {'rms_s_err':>10}  "
          f"{'refits':>6}  {'sec':>5}  {'final_s':>7}  {'score':>7}")
    best = None
    for W in W_grid:
        sdm, sdd = floors[W]
        for k in k_grid:
            res = run_tracker(g["log_pdf"], g["sampler"], g["theta_seed"],
                              df, X_drift, s_truth, W, k, sdm, sdd, sim_mult)
            score = res["rms_s_error"] + lambda_cost * res["refits"] / n_run * 1000
            rows.append((W, k, res, score))
            star = ""
            if best is None or score < best[-1]:
                best = (W, k, res, score)
                star = "  <-"
            print(f"  {W:>6}  {k:>7.1f}  {res['rms_s_error']:>10.4f}  "
                  f"{res['refits']:>6d}  {res['elapsed']:>5.1f}  "
                  f"{res['s_path_final']:>7.3f}  {score:>7.4f}{star}")

    W_best, k_best, res_best, _ = best
    config = {
        "df": df,
        "W": W_best,
        "alpha": 1.0 / W_best,
        "k_sigma": k_best,
        "sim_mult": sim_mult,
        "M": sim_mult * W_best,
        "refit_floor": W_best,
        "burn_in": max(1, W_best // 20),
        "sd_mean_empirical": floors[W_best][0],
        "sd_d2_empirical": floors[W_best][1],
        "theta_seed": g["theta_seed"].tolist(),
        "rms_L_pool": g["rms_L_pool"],
        "drift_score": res_best["rms_s_error"],
    }
    out = Path(__file__).resolve().parent / "luck_config.json"
    out.write_text(json.dumps(config, indent=2))
    print(f"\nselected: W={W_best}, k_sigma={k_best}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
