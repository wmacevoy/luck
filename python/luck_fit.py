"""Generic luck-based parameter fitting.

Implements the numerical-luck framework from computation.tex / numluck.sci /
numlucksetup.sci as a parameter-fitting loop.  For a parameterized family
p(.|theta) supplied as (log_pdf, sampler), the chapter defines

    L_theta(x) = Pr_{X ~ p(.|theta)} ( log p(X|theta) >= log p(x|theta) )

and estimates it from a simulated set S = {x_1, ..., x_M} ~ p(.|theta) by

    ell_theta(x) = (#S more probable than x)/M + (1/2)(#S equally probable)/M

with E[ell] = L and Var(ell) <= L(1-L)/M <= 1/(4M).  Under a correct fit
the random variable L_theta(X_data) is uniform on [0,1].

Fit objective = RMS deviation of the sorted empirical {ell_theta(x_i)} from
the uniform quantiles (i - 1/2)/n_data (Cramer-von Mises distance).

Common random numbers (fixed simulation seed across optimizer steps)
keep the objective smooth in theta so vanilla Nelder-Mead converges.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln


class LuckFit:
    def __init__(self, log_pdf, sampler, theta0, n_sim: int = 4000,
                 sim_seed: int = 0):
        self.log_pdf = log_pdf
        self.sampler = sampler
        self.theta = np.asarray(theta0, dtype=float).copy()
        self.n_sim = int(n_sim)
        self.sim_seed = int(sim_seed)
        self.X: list[np.ndarray] = []
        self.fit_result = None

    def add_samples(self, X) -> "LuckFit":
        for x in np.atleast_2d(np.asarray(X, dtype=float)):
            self.X.append(x.copy())
        return self

    def ell_pair(self, theta, X_data):
        """Return (ell, 1 - ell), each computed from rank counts so neither
        extreme loses precision via 1 - small.  ell is the chapter's
        count-based estimator of L = |Omega| + 1/2 |omega|."""
        rng = np.random.default_rng(self.sim_seed)
        X_sim = self.sampler(theta, self.n_sim, rng)
        nlp_sim = -self.log_pdf(theta, X_sim)
        nlp_data = -self.log_pdf(theta, X_data)
        nlp_sim_sorted = np.sort(nlp_sim)
        M = nlp_sim_sorted.size
        idx_lo = np.searchsorted(nlp_sim_sorted, nlp_data, side="left")
        idx_hi = np.searchsorted(nlp_sim_sorted, nlp_data, side="right")
        num_ell = idx_lo + 0.5 * (idx_hi - idx_lo)            # |A| + 1/2 |a|
        num_comp = (M - idx_hi) + 0.5 * (idx_hi - idx_lo)     # complement
        return num_ell / M, num_comp / M

    def ell(self, theta, X_data) -> np.ndarray:
        e, _ = self.ell_pair(theta, X_data)
        return e

    def coin_luck(self, theta, X_data) -> np.ndarray:
        """C = log2(ell / (1 - ell)).  Coin-luck: how many coin-flips ahead
        of (C>0, rare/improbable tail) or behind (C<0, over-common tail)
        the data is, relative to typical.  Floored at |C| <= log2(2M) --
        the physical resolution given M simulated samples (one half-sample)."""
        e, comp = self.ell_pair(theta, X_data)
        floor = 0.5 / self.n_sim
        return np.log2(np.maximum(e, floor)) - np.log2(np.maximum(comp, floor))

    def rms_error(self, theta) -> float:
        """Midrank-corrected CvM-style distance from uniform.

        Standard target (i + 1/2)/n assumes distinct ell values; under a
        misspecified model, ell saturates at 0 or 1 with multiplicity, so
        each tied block is collapsed to a single midrank target equal to
        the mean of the would-be targets within the block.  Without this,
        ties contribute irreducible variance that the optimizer chases
        but cannot reduce.
        """
        X = np.asarray(self.X)
        ell_sorted = np.sort(self.ell(theta, X))
        n = ell_sorted.size
        target = np.empty(n)
        i = 0
        while i < n:
            j = i
            while j + 1 < n and ell_sorted[j + 1] == ell_sorted[i]:
                j += 1
            target[i:j + 1] = (i + j + 1) / (2 * n)
            i = j + 1
        return float(np.sqrt(np.mean((ell_sorted - target) ** 2)))

    def fit(self, method: str = "Nelder-Mead", **kwargs) -> "LuckFit":
        kwargs.setdefault("options",
                          {"xatol": 1e-4, "fatol": 1e-6, "maxiter": 4000})
        self.fit_result = minimize(self.rms_error, self.theta,
                                   method=method, **kwargs)
        self.theta = np.asarray(self.fit_result.x, dtype=float)
        return self


def mvt_family(mu_shape, Sigma_shape):
    """Build (log_pdf, sampler) for a multivariate-t with theta = (s, log_nu).

    mean = mu_shape, scale matrix = s^2 * Sigma_shape, df = exp(log_nu).
    """
    mu_shape = np.asarray(mu_shape, dtype=float)
    k = mu_shape.size
    L_shape = np.linalg.cholesky(Sigma_shape)
    log_det_shape = 2.0 * np.sum(np.log(np.diag(L_shape)))

    def log_pdf(theta, X):
        s, log_nu = theta
        nu = float(np.exp(log_nu))
        Xc = X - mu_shape
        Z = np.linalg.solve(L_shape, Xc.T) / s
        r2 = np.einsum("ij,ij->j", Z, Z)
        log_const = (gammaln((nu + k) / 2) - gammaln(nu / 2)
                     - 0.5 * k * np.log(np.pi * nu))
        log_det = log_det_shape + 2.0 * k * np.log(s)
        return log_const - 0.5 * log_det \
               - 0.5 * (nu + k) * np.log1p(r2 / nu)

    def sampler(theta, n, rng):
        s, log_nu = theta
        nu = float(np.exp(log_nu))
        g = rng.chisquare(nu, size=(n, 1)) / nu
        Z = rng.standard_normal((n, k))
        return mu_shape + s * (Z @ L_shape.T) / np.sqrt(g)

    return log_pdf, sampler


def _gauss(rng, mu, L, n):
    return mu + (L @ rng.standard_normal((mu.size, n))).T


def _t(rng, mu, L, nu, n):
    g = rng.chisquare(nu, size=(n, 1)) / nu
    return mu + (L @ (rng.standard_normal((mu.size, n)) / np.sqrt(g.T))).T


if __name__ == "__main__":
    rng = np.random.default_rng(20260430)
    k = 12
    A = rng.standard_normal((k, k))
    Sigma_true = A @ A.T / k + 0.5 * np.eye(k)
    mu_true = rng.standard_normal(k)
    L_true = np.linalg.cholesky(Sigma_true)

    n = 1500

    def run(name, X):
        mu_hat = X.mean(axis=0)
        Sigma_hat = np.cov(X.T, bias=False)
        log_pdf, sampler = mvt_family(mu_hat, Sigma_hat)
        lf = LuckFit(log_pdf, sampler, theta0=[1.0, np.log(60)], n_sim=4000)
        lf.add_samples(X)
        rms_init = lf.rms_error(lf.theta)
        lf.fit()
        s_hat, log_nu_hat = lf.theta
        rms_final = lf.rms_error(lf.theta)
        print(f"\n=== {name} (n={len(X)}, k={k}) ===")
        print(f"  init theta=(s=1.000, nu=60.0)   rms_L = {rms_init:.5f}")
        print(f"  fit  theta=(s={s_hat:.3f},   nu={np.exp(log_nu_hat):.2f})  "
              f"rms_L = {rms_final:.5f}")
        print(f"  optimizer: {lf.fit_result.message}  "
              f"nfev={lf.fit_result.nfev}")
        return lf

    run("Gaussian data", _gauss(rng, mu_true, L_true, n))
    run("t(nu=5) data",  _t(rng, mu_true, L_true, 5.0, n))
    run("t(nu=3) data",  _t(rng, mu_true, L_true, 3.0, n))

    print("\n=== sanity: Gaussian-only family on t(5) data ===")
    nu_true = 5.0
    X_t = _t(rng, mu_true, L_true, nu_true, n)
    mu_hat = X_t.mean(axis=0)
    Sigma_hat = np.cov(X_t.T, bias=False)
    log_pdf, sampler = mvt_family(mu_hat, Sigma_hat)

    def gauss_log_pdf(theta, X):
        return log_pdf((theta[0], np.log(1e6)), X)

    def gauss_sampler(theta, n_, rng_):
        return sampler((theta[0], np.log(1e6)), n_, rng_)

    lf = LuckFit(gauss_log_pdf, gauss_sampler, theta0=[1.0], n_sim=4000)
    lf.add_samples(X_t)
    lf.fit()
    print(f"  Gaussian-restricted fit: s={lf.theta[0]:.3f}  "
          f"rms_L = {lf.rms_error(lf.theta):.5f}  "
          f"(big rms = model class is wrong, not a parameter problem)")
