"""Online luck tracker: EWMA monitor + event-triggered local refit.

Per-step work is cheap: one log_pdf evaluation, one searchsorted lookup,
two EWMA updates.  When the windowed L statistics drift outside the noise
floor of a correctly-fit model, a short local Nelder-Mead refit runs on
the current data window starting from the present theta.

Window size W = 24 * df^2 is chosen so the per-sample numerical-luck
noise sd <= 1/(2 sqrt(W)) is O(1/df) -- comparable to one CDF "tick"
under the model -- and so a Gaussian shape (df(df+1)/2 parameters) is
well-determined within W.

Trigger rule (under a correct fit, ell ~ Uniform(0,1)):
  E[ell]            = 1/2          sd of EWMA mean ~ sqrt(1/(12 W))
  E[(ell-1/2)^2]    = 1/12         sd of EWMA-var-proxy ~ sqrt(c/W),  c ~ 1/180
A refit fires when either EWMA leaves a k_sigma-wide band, or after a
hard floor of W / refit_floor_div steps without one (slow-drift safety
net).  Refits refresh the simulation buffer at the new theta and reset
the EWMAs to their target values to avoid spurious re-triggering during
burn-in.
"""

from __future__ import annotations

from collections import deque

import numpy as np
from scipy.optimize import minimize


class LuckTrackerOnline:
    def __init__(self, log_pdf, sampler, theta0, df: int,
                 k_sigma: float = 3.0,
                 sim_mult: int = 4,
                 refit_floor_div: int = 1,
                 max_refit_iters: int = 30,
                 burn_in_div: int = 20,
                 sim_seed: int = 0):
        self.log_pdf = log_pdf
        self.sampler = sampler
        self.theta = np.asarray(theta0, dtype=float).copy()
        self.df = int(df)
        self.W = 24 * self.df ** 2
        self.alpha = 1.0 / self.W
        self.k_sigma = float(k_sigma)
        self.refit_floor = max(1, self.W // refit_floor_div)
        self.max_refit_iters = int(max_refit_iters)
        self.burn_in = max(1, self.W // burn_in_div)
        self.M = sim_mult * self.W
        self._rng = np.random.default_rng(sim_seed)

        self.t = 0
        self.window: deque[np.ndarray] = deque(maxlen=self.W)
        self.steps_since_refit = 0
        self.ewma_mean = 0.5
        self.ewma_d2 = 1.0 / 12.0
        self.refits = 0
        self.last_refit_t = 0
        self.history: list[dict] = []

        self._sim_nlp_sorted: np.ndarray | None = None
        self._refresh_sim_buffer()

    def _refresh_sim_buffer(self):
        X_sim = self.sampler(self.theta, self.M, self._rng)
        nlp = -self.log_pdf(self.theta, X_sim)
        self._sim_nlp_sorted = np.sort(nlp)

    def _ell_pair(self, x: np.ndarray) -> tuple[float, float]:
        """Return (L, U) for sample x against the current sim buffer.
        Both come from rank counts; neither is computed via 1 - other,
        so precision is preserved at both extremes (and no clipping)."""
        nlp_x = float(-self.log_pdf(self.theta, x.reshape(1, -1))[0])
        a = self._sim_nlp_sorted
        M = a.size
        idx_lo = int(np.searchsorted(a, nlp_x, side="left"))
        idx_hi = int(np.searchsorted(a, nlp_x, side="right"))
        L = (idx_lo + 0.5 * (idx_hi - idx_lo)) / M
        U = ((M - idx_hi) + 0.5 * (idx_hi - idx_lo)) / M
        return L, U

    def _ell(self, x: np.ndarray) -> float:
        L, _ = self._ell_pair(x)
        return L

    def update(self, x) -> float:
        x = np.asarray(x, dtype=float).reshape(self.df)
        self.t += 1
        self.window.append(x.copy())
        ell = self._ell(x)
        d = ell - 0.5
        a = self.alpha
        self.ewma_mean = (1 - a) * self.ewma_mean + a * ell
        self.ewma_d2 = (1 - a) * self.ewma_d2 + a * d * d
        self.steps_since_refit += 1

        triggered_by = None
        if self.t - self.last_refit_t > self.burn_in:
            sd_mean = np.sqrt((1.0 / 12.0) / (2.0 * self.W))
            sd_d2 = np.sqrt((1.0 / 180.0) / (2.0 * self.W))
            if abs(self.ewma_mean - 0.5) > self.k_sigma * sd_mean:
                triggered_by = "mean"
            elif abs(self.ewma_d2 - 1.0 / 12.0) > self.k_sigma * sd_d2:
                triggered_by = "var"

        if triggered_by is not None or self.steps_since_refit >= self.refit_floor:
            self._local_refit(triggered_by or "floor")

        return ell

    def _local_refit(self, reason: str):
        if len(self.window) < max(self.df + 2, 30):
            self.steps_since_refit = 0
            return

        X_data = np.asarray(self.window)
        self._refresh_sim_buffer()
        sim_seed_for_loss = int(self.t)

        def loss(theta):
            rng = np.random.default_rng(sim_seed_for_loss)
            X_sim = self.sampler(theta, self.M, rng)
            nlp_sim = np.sort(-self.log_pdf(theta, X_sim))
            nlp_data = -self.log_pdf(theta, X_data)
            idx_lo = np.searchsorted(nlp_sim, nlp_data, side="left")
            idx_hi = np.searchsorted(nlp_sim, nlp_data, side="right")
            ell = np.sort(0.5 * (idx_lo + idx_hi) / nlp_sim.size)
            n = ell.size
            target = np.empty(n)
            i = 0
            while i < n:
                j = i
                while j + 1 < n and ell[j + 1] == ell[i]:
                    j += 1
                target[i:j + 1] = (i + j + 1) / (2 * n)
                i = j + 1
            return float(np.sqrt(np.mean((ell - target) ** 2)))

        scale = np.maximum(np.abs(self.theta), 1e-3) * 0.05
        simplex = np.vstack([self.theta] + [self.theta + scale * e
                                            for e in np.eye(self.theta.size)])
        res = minimize(loss, self.theta, method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-4,
                                "maxiter": self.max_refit_iters,
                                "initial_simplex": simplex})
        self.theta = np.asarray(res.x, dtype=float)
        self.refits += 1
        self.last_refit_t = self.t
        self.steps_since_refit = 0
        self.ewma_mean = 0.5
        self.ewma_d2 = 1.0 / 12.0
        self.history.append({"t": self.t, "reason": reason,
                             "theta": self.theta.copy(),
                             "rms_L": res.fun, "nfev": res.nfev})
        self._refresh_sim_buffer()


if __name__ == "__main__":
    from luck_fit import mvt_family

    rng = np.random.default_rng(20260430)
    k = 12

    A = rng.standard_normal((k, k))
    Sigma_base = A @ A.T / k + 0.5 * np.eye(k)
    mu_true = rng.standard_normal(k)
    L_base = np.linalg.cholesky(Sigma_base)

    n_burn = 1500
    X_burn = mu_true + (L_base @ rng.standard_normal((k, n_burn))).T
    mu_hat = X_burn.mean(axis=0)
    Sigma_hat = np.cov(X_burn.T, bias=False)

    log_pdf, sampler = mvt_family(mu_hat, Sigma_hat)

    n_drift = 8000

    def s_true_at(t):
        if t < 2000:
            return 1.0
        if t < 5000:
            return 1.0 + 0.5 * (t - 2000) / 3000
        return 1.5

    tracker = LuckTrackerOnline(log_pdf, sampler,
                                theta0=[1.0, np.log(60)], df=k,
                                k_sigma=3.0)

    sd_mean = np.sqrt((1/12)/(2*tracker.W))
    sd_d2   = np.sqrt((1/180)/(2*tracker.W))
    print(f"W = 24 * df^2 = {tracker.W}")
    print(f"refit floor every {tracker.refit_floor} steps "
          f"(burn-in {tracker.burn_in})")
    print(f"sim buffer M = {tracker.M}")
    print(f"3-sigma trigger bands: |mean-1/2| > {3*sd_mean:.4f},  "
          f"|d^2-1/12| > {3*sd_d2:.4f}\n")

    print(f"  t    s_true   s_hat   nu_hat      ewma(mean) ewma(d^2)  refits  last")
    print(f"  {'-'*82}")
    for t in range(1, n_drift + 1):
        s_t = s_true_at(t)
        x = mu_true + s_t * (L_base @ rng.standard_normal(k))
        tracker.update(x)
        if t % 500 == 0:
            s_hat, log_nu_hat = tracker.theta
            last = tracker.history[-1]["reason"] if tracker.history else "-"
            print(f"  {t:5d}  {s_t:.3f}   {s_hat:.3f}   "
                  f"{np.exp(log_nu_hat):8.1f}    {tracker.ewma_mean:.4f}    "
                  f"{tracker.ewma_d2:.4f}     {tracker.refits:3d}    {last}")

    print(f"\ntotal refits: {tracker.refits} over {n_drift} samples "
          f"(amortized 1 per {n_drift // max(1, tracker.refits)} samples)")
    refit_reasons = [h["reason"] for h in tracker.history]
    from collections import Counter
    print(f"refit reasons: {dict(Counter(refit_reasons))}")
    print(f"final s_hat = {tracker.theta[0]:.3f}  "
          f"(s_true at end = {s_true_at(n_drift):.3f})")
