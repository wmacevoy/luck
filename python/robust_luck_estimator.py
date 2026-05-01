"""Robust luck estimator: multivariate-t MLE via EM, luck against t-shell.

For elliptical-but-heavier-than-Gaussian data, a single rescaling of the
sample covariance cannot fix the radial distribution -- the *shape* of the
shell is wrong, not just its size.  Fitting a multivariate t with
(mu, Sigma, nu) gets both right:

  R^2 = (x - mu)^T Sigma^{-1} (x - mu) ~ k * F(k, nu)

and the luck shell radius is

  c_t(k, nu) = E[R] = Gamma((k+1)/2)/Gamma(k/2)
                      * sqrt(nu) * Gamma((nu-1)/2)/Gamma(nu/2)

which tends to sqrt(k - 1/2) as nu -> infinity, recovering the Gaussian
luck shell.

EM update (per refit):
  E:  w_i = (nu + k) / (nu + r_i^2)
  M:  mu     = sum w_i x_i / sum w_i
      Sigma  = (1/n) sum w_i (x_i - mu)(x_i - mu)^T
      nu     : root of digamma equation, by bisection
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln, digamma
from scipy.optimize import brentq


class RobustLuckEstimator:
    def __init__(self, dim: int, em_iters: int = 30, ridge: float = 1e-9,
                 nu_min: float = 2.5, nu_max: float = 500.0):
        self.k = int(dim)
        self.em_iters = int(em_iters)
        self.ridge = float(ridge)
        self.nu_min = float(nu_min)
        self.nu_max = float(nu_max)
        self.n = 0
        self._X: list[np.ndarray] = []
        self.mu = np.zeros(self.k)
        self.Sigma = np.eye(self.k)
        self.nu = self.nu_max

    def update(self, x) -> "RobustLuckEstimator":
        x = np.asarray(x, dtype=float).reshape(self.k)
        self._X.append(x.copy())
        self.n += 1
        if self.n >= max(2 * self.k, self.k + 5):
            self._refit()
        return self

    def _r2(self, X, mu, Sigma):
        L = np.linalg.cholesky(Sigma + self.ridge * np.eye(self.k))
        Z = np.linalg.solve(L, (X - mu).T)
        return np.einsum("ij,ij->j", Z, Z)

    def _refit(self):
        X = np.asarray(self._X)
        n, k = X.shape
        mu = X.mean(axis=0)
        Sigma = np.cov(X.T, bias=False) + self.ridge * np.eye(k)
        nu = self.nu

        for _ in range(self.em_iters):
            r2 = self._r2(X, mu, Sigma)
            w = (nu + k) / (nu + r2)

            wsum = w.sum()
            mu_new = (w[:, None] * X).sum(axis=0) / wsum
            Xc = X - mu_new
            Sigma_new = (Xc.T * w) @ Xc / n + self.ridge * np.eye(k)

            r2_new = self._r2(X, mu_new, Sigma_new)
            w_new = (nu + k) / (nu + r2_new)
            const = 1.0 + np.mean(np.log(w_new) - w_new) \
                    + digamma((nu + k) / 2) - np.log((nu + k) / 2)

            def g(nu_):
                return -digamma(nu_ / 2) + np.log(nu_ / 2) + const

            try:
                nu_new = brentq(g, self.nu_min, self.nu_max,
                                xtol=1e-4, maxiter=100)
            except ValueError:
                nu_new = self.nu_max if g(self.nu_max) > 0 else self.nu_min

            if (np.linalg.norm(mu_new - mu) < 1e-7 and
                    np.linalg.norm(Sigma_new - Sigma) < 1e-7 and
                    abs(nu_new - nu) < 1e-3):
                mu, Sigma, nu = mu_new, Sigma_new, nu_new
                break
            mu, Sigma, nu = mu_new, Sigma_new, nu_new

        self.mu, self.Sigma, self.nu = mu, Sigma, nu

    def shell(self) -> float:
        k, nu = self.k, self.nu
        if nu > 200:
            return float(np.sqrt(k - 0.5))
        log_c = (gammaln((k + 1) / 2) - gammaln(k / 2)
                 + 0.5 * np.log(nu)
                 + gammaln((nu - 1) / 2) - gammaln(nu / 2))
        return float(np.exp(log_c))

    def z_luck(self, x) -> np.ndarray:
        """z_L = r - shell(nu, k), the t-distribution analogue of the
        Gaussian z-luck from normal.tex."""
        X = np.atleast_2d(np.asarray(x, dtype=float))
        if X.shape[1] != self.k and X.shape[0] == self.k:
            X = X.T
        r = np.sqrt(self._r2(X, self.mu, self.Sigma))
        return r - self.shell()


def _gauss_samples(rng, mu, L, n):
    return mu + (L @ rng.standard_normal((mu.size, n))).T


def _t_samples(rng, mu, L, nu, n):
    g = rng.chisquare(nu, size=(n, 1)) / nu
    return mu + (L @ (rng.standard_normal((mu.size, n)) / np.sqrt(g.T))).T


if __name__ == "__main__":
    from luck_estimator import LuckEstimator

    rng = np.random.default_rng(20260430)
    k = 12
    A = rng.standard_normal((k, k))
    Sigma_true = A @ A.T / k + 0.5 * np.eye(k)
    mu_true = rng.standard_normal(k)
    L_true = np.linalg.cholesky(Sigma_true)

    n_train, n_test = 2000, 20000

    def run(name, X_train, X_test, true_shell):
        gauss_est = LuckEstimator(k)
        for x in X_train:
            gauss_est.update(x)
        rob_est = RobustLuckEstimator(k)
        for x in X_train:
            rob_est.update(x)

        z_gauss = gauss_est.z_luck(X_test)
        z_rob = rob_est.z_luck(X_test)

        def stats(z):
            return f"mean={z.mean():+.4f}  std={z.std():.4f}  rms={np.sqrt((z**2).mean()):.4f}"

        print(f"\n=== {name} ===  (true t-shell = {true_shell:.4f})")
        print(f"  Gaussian luck-fit:  {stats(z_gauss)}")
        print(f"     -> shell used = {gauss_est.shell:.4f}")
        print(f"  t-MLE luck-fit:     {stats(z_rob)}")
        print(f"     -> nu_hat={rob_est.nu:.2f}  shell used={rob_est.shell():.4f}")

    X_tr = _gauss_samples(rng, mu_true, L_true, n_train)
    X_te = _gauss_samples(rng, mu_true, L_true, n_test)
    run("Gaussian data", X_tr, X_te,
        true_shell=float(np.sqrt(k - 0.5)))

    nu_true = 5.0
    X_tr = _t_samples(rng, mu_true, L_true, nu_true, n_train)
    X_te = _t_samples(rng, mu_true, L_true, nu_true, n_test)
    log_c = (gammaln((k + 1) / 2) - gammaln(k / 2)
             + 0.5 * np.log(nu_true)
             + gammaln((nu_true - 1) / 2) - gammaln(nu_true / 2))
    run(f"t(nu={nu_true}) data", X_tr, X_te,
        true_shell=float(np.exp(log_c)))

    nu_true = 3.5
    X_tr = _t_samples(rng, mu_true, L_true, nu_true, n_train)
    X_te = _t_samples(rng, mu_true, L_true, nu_true, n_test)
    log_c = (gammaln((k + 1) / 2) - gammaln(k / 2)
             + 0.5 * np.log(nu_true)
             + gammaln((nu_true - 1) / 2) - gammaln(nu_true / 2))
    run(f"t(nu={nu_true}) data", X_tr, X_te,
        true_shell=float(np.exp(log_c)))
