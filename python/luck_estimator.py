"""Online multivariate-normal estimator that minimizes squared luck error.

Standard sample-covariance estimators target the *center* of the predictive
distribution.  In k dimensions, samples actually concentrate on a shell of
radius c = sqrt(k - 1/2).  This estimator picks (mu, Sigma) so that the
empirical luck statistic z = r - c is as close to zero as possible in the
mean-squared sense, where r = sqrt((x - mu)^T Sigma^{-1} (x - mu)).

Update rule:
  - mu      : running sample mean (Welford).
  - shape S : running sample covariance (Welford M2 / (n-1)).
  - scale s : minimizes (1/n) sum_i (u_i / s - c)^2 with u_i computed from S,
              giving the closed form  s = mean(u^2) / (c * mean(u)).
  - Sigma   : s^2 * S.

For a true Gaussian, mean(u) -> c and mean(u^2) -> k, so s -> k/(k - 1/2)
which is ~1 for moderate k -- the luck rescaling is a small correction on
clean data and a meaningful one on heavy-tailed or contaminated data.
"""

from __future__ import annotations

import numpy as np


class LuckEstimator:
    def __init__(self, dim: int, ridge: float = 1e-9):
        self.k = int(dim)
        self.shell = float(np.sqrt(self.k - 0.5))
        self.n = 0
        self.mu = np.zeros(self.k)
        self._M2 = np.zeros((self.k, self.k))
        self._X: list[np.ndarray] = []
        self.scale = 1.0
        self.ridge = float(ridge)

    def update(self, x) -> "LuckEstimator":
        x = np.asarray(x, dtype=float).reshape(self.k)
        self.n += 1
        d1 = x - self.mu
        self.mu = self.mu + d1 / self.n
        d2 = x - self.mu
        self._M2 = self._M2 + np.outer(d1, d2)
        self._X.append(x.copy())
        if self.n >= max(2, self.k + 1):
            self._refit_scale()
        return self

    def _shape(self) -> np.ndarray:
        if self.n < 2:
            return np.eye(self.k)
        return self._M2 / (self.n - 1) + self.ridge * np.eye(self.k)

    def _refit_scale(self) -> None:
        S = self._shape()
        L = np.linalg.cholesky(S)
        X = np.asarray(self._X) - self.mu
        Z = np.linalg.solve(L, X.T)
        u = np.sqrt(np.einsum("ij,ij->j", Z, Z))
        um = u.mean()
        if um > 0:
            self.scale = float((u * u).mean() / (self.shell * um))

    @property
    def Sigma(self) -> np.ndarray:
        return (self.scale ** 2) * self._shape()

    def z_luck(self, x) -> np.ndarray:
        """Return z_L = r - sqrt(k - 1/2) for one or many samples
        (chapter normal.tex, def. of z-luck z_L)."""
        X = np.atleast_2d(np.asarray(x, dtype=float))
        if X.shape[1] != self.k and X.shape[0] == self.k:
            X = X.T
        L = np.linalg.cholesky(self.Sigma)
        Z = np.linalg.solve(L, (X - self.mu).T)
        r = np.sqrt(np.einsum("ij,ij->j", Z, Z))
        return r - self.shell


def _sample_gaussian(rng, mu, L, n):
    return mu + (L @ rng.standard_normal((mu.size, n))).T


def _z_luck_stats(mu_hat, Sigma_hat, X, shell):
    L = np.linalg.cholesky(Sigma_hat)
    Z = np.linalg.solve(L, (X - mu_hat).T)
    r = np.sqrt(np.einsum("ij,ij->j", Z, Z))
    return r - shell


if __name__ == "__main__":
    rng = np.random.default_rng(20260430)

    k = 12
    A = rng.standard_normal((k, k))
    Sigma_true = A @ A.T / k + 0.5 * np.eye(k)
    mu_true = rng.standard_normal(k)
    L_true = np.linalg.cholesky(Sigma_true)

    n_train = 2000
    n_test = 20000

    X_train = _sample_gaussian(rng, mu_true, L_true, n_train)
    X_test = _sample_gaussian(rng, mu_true, L_true, n_test)

    est = LuckEstimator(k)
    for x in X_train:
        est.update(x)

    z_fit = est.z_luck(X_test)

    mu_naive = X_train.mean(axis=0)
    Sigma_naive = np.cov(X_train.T, bias=False)
    z_naive = _z_luck_stats(mu_naive, Sigma_naive, X_test, est.shell)

    z_oracle = _z_luck_stats(mu_true, Sigma_true, X_test, est.shell)

    target_std = float(np.sqrt(0.5))

    def report(name, z):
        print(f"  {name:12s}  mean={z.mean():+.4f}  "
              f"std={z.std():.4f}  rms={np.sqrt((z**2).mean()):.4f}")

    print(f"k = {k}    shell c = sqrt(k - 1/2) = {est.shell:.4f}")
    print(f"target (Gaussian, large n):  mean=+0.0000  std~{target_std:.4f}\n")
    print(f"out-of-sample z-luck (n_test={n_test}):")
    report("luck-fit",   z_fit)
    report("sample-cov", z_naive)
    report("oracle",     z_oracle)
    print(f"\nfitted scale s = {est.scale:.5f}  "
          f"(asymptotic target for k={k}: {k / (k - 0.5):.5f})")

    rng2 = np.random.default_rng(7)
    nu = 5
    g = rng2.chisquare(nu, size=(n_train, 1)) / nu
    X_heavy_train = mu_true + (L_true @ (rng2.standard_normal((k, n_train)) /
                                          np.sqrt(g.T))).T
    g2 = rng2.chisquare(nu, size=(n_test, 1)) / nu
    X_heavy_test = mu_true + (L_true @ (rng2.standard_normal((k, n_test)) /
                                          np.sqrt(g2.T))).T

    est_h = LuckEstimator(k)
    for x in X_heavy_train:
        est_h.update(x)

    print(f"\nheavy-tailed (multivariate t, nu={nu}) -- model is misspecified:")
    report("luck-fit",   est_h.z_luck(X_heavy_test))
    mu_h = X_heavy_train.mean(axis=0)
    S_h = np.cov(X_heavy_train.T, bias=False)
    report("sample-cov", _z_luck_stats(mu_h, S_h, X_heavy_test, est_h.shell))
    print(f"fitted scale s = {est_h.scale:.5f}  "
          f"(>1 means luck-fit inflated Sigma vs sample cov)")
