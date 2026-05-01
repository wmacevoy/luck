"""Verify L, U, C carry signal at the extremes (no silent clipping)."""

import numpy as np

from luck_fit import LuckFit, mvt_family


def main():
    rng = np.random.default_rng(0)
    k = 4
    Sigma = np.eye(k) * 4.0
    mu = np.zeros(k)
    L_chol = np.linalg.cholesky(Sigma)

    n = 500
    X_data = mu + (rng.standard_normal((n, k)) @ L_chol.T)

    Sigma_tight = 0.25 * np.eye(k)
    log_pdf, sampler = mvt_family(mu, Sigma_tight)
    lf = LuckFit(log_pdf, sampler, theta0=[1.0, np.log(60)], n_sim=2000)

    ell, comp = lf.ell_pair(lf.theta, X_data)
    C = lf.coin_luck(lf.theta, X_data)

    print("Misspecified model (data sd=2 vs model sd=0.5), n=500")
    print(f"  ell:    min={ell.min():.4g}  max={ell.max():.4g}  "
          f"frac at exact 1.0 = {(ell == 1.0).mean():.2%}")
    print(f"  1-ell:  min={comp.min():.4g}  max={comp.max():.4g}  "
          f"frac at exact 0.0 = {(comp == 0.0).mean():.2%}")
    print(f"  C: min={C.min():+.3f}  median={np.median(C):+.3f}  max={C.max():+.3f}")
    print(f"     log2(2M) cap = {np.log2(2 * lf.n_sim):.3f}")
    print(f"  ell+(1-ell) identity (should be 1.0): max |sum-1| = {np.abs(ell+comp-1).max():.2e}")

    Sigma_correct = Sigma
    log_pdf2, sampler2 = mvt_family(mu, Sigma_correct)
    lf2 = LuckFit(log_pdf2, sampler2, theta0=[1.0, np.log(60)], n_sim=2000)
    ell2, _ = lf2.ell_pair(lf2.theta, X_data)
    C2 = lf2.coin_luck(lf2.theta, X_data)
    print("\nCorrectly specified model, same data:")
    print(f"  ell: min={ell2.min():.4g}  max={ell2.max():.4g}")
    print(f"  C:   min={C2.min():+.3f}  median={np.median(C2):+.3f}  "
          f"max={C2.max():+.3f}")
    print(f"  mean|C| = {np.abs(C2).mean():.3f}  "
          f"(vs misspecified mean|C| = {np.abs(C).mean():.3f})")


if __name__ == "__main__":
    main()
