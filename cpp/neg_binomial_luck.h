#pragma once

#include "luck.h"
#include <cmath>
#include <limits>

namespace luck {

  // ---------------------------------------------------------------------------
  // Negative binomial: number of failures k = 0, 1, 2, ... before the r-th
  // success in a sequence of Bernoulli(p) trials. r is the target success
  // count (integer >= 1), q = 1 - p.
  //   pmf(k) = C(k+r-1, k) * p^r * q^k
  // r = 1 is the geometric distribution.
  // ---------------------------------------------------------------------------
  template <typename R>
  R neg_binomial_log_pmf(int r, R p, R q, R x) {
    return std::lgamma(x + R(r))
         - std::lgamma(R(r))
         - std::lgamma(x + R(1))
         + R(r) * std::log(p)
         + x    * std::log(q);
  }

  template <typename R>
  R neg_binomial_pmf(int r, R p, R q, int k) {
    return std::exp(neg_binomial_log_pmf<R>(r, p, q, R(k)));
  }

  // ---------------------------------------------------------------------------
  // Negative-binomial luck mass at observation k of NB(r, p):
  //   less  = sum_{j: pmf(j) <  pmf(k)} pmf(j)
  //   equal = sum_{j: pmf(j) == pmf(k)} pmf(j)   (= pmf(k), or 2*pmf(k) at a
  //                                               tied integer-mode peak)
  //   more  = sum_{j: pmf(j) >  pmf(k)} pmf(j)
  //
  // Same band-walk pattern as poisson_mass: the pmf is unimodal and the ratio
  //   pmf(j+1)/pmf(j) = q * (j + r) / (j + 1)
  // is decreasing in j, so we walk from k across the mode -- no root solve --
  // accumulating `more` term by term, and stop when the pmf falls back below
  // pmf(k). `less` = 1 - equal - more, so the infinite tail is never summed.
  // ---------------------------------------------------------------------------
  template <typename R>
  mass<R> neg_binomial_mass(int r, R p, R q, int k) {
    if (k < 0 || r < 1) return { R(0), R(0), R(0) };
    if (p <= R(0))      return { R(0), R(0), R(0) };          // never succeeds
    if (p >= R(1)) {                                          // pmf(0)=1, rest 0
      if (k == 0)       return { R(0), R(1), R(0) };
      return { R(0), R(0), R(1) };
    }

    // Mode of the unimodal pmf. The ratio q*(k+r)/(k+1) equals 1 at the real
    // value (q*r - 1)/p; ceiling that gives the peak (or 0 if it's negative),
    // and an exact integer value gives a tied pair at (mode-1, mode).
    R mode_real = (q * R(r) - R(1)) / p;
    int  mode      = (mode_real < R(0)) ? 0 : int(std::floor(mode_real)) + 1;
    bool mode_pair = (mode >= 1) && (R(mode - 1) == mode_real);

    R pmf_k = neg_binomial_pmf<R>(r, p, q, k);

    if (k == mode || (mode_pair && k == mode - 1)) {
      R eq = (mode_pair ? R(2) : R(1)) * pmf_k;
      return { R(1) - eq, eq, R(0) };
    }

    R eps    = R(64) * std::numeric_limits<R>::epsilon();
    R cutoff = pmf_k - eps * pmf_k;
    R more   = R(0);
    R pmf_el = R(0);

    if (k < mode) {
      // rising flank: walk up, pmf(j+1) = pmf(j) * q*(j+r)/(j+1).
      R pj = pmf_k;
      for (int j = k; ; ++j) {
        pj *= q * R(j + r) / R(j + 1);     // pj = pmf(j+1)
        if (pj < cutoff) break;
        pmf_el = pj;
        more  += pj;
      }
    } else {
      // falling flank: walk down, pmf(j-1) = pmf(j) * j / (q*(j+r-1)).
      R pj = pmf_k;
      for (int j = k; j > 0; --j) {
        pj *= R(j) / (q * R(j + r - 1));   // pj = pmf(j-1)
        if (pj < cutoff) break;
        pmf_el = pj;
        more  += pj;
      }
    }

    // el is the far band endpoint; tie -> equal, otherwise it stays in more.
    int same = std::abs(pmf_el - pmf_k) <= eps * pmf_k;
    if (same) more -= pmf_el;

    R equal = R(1 + same) * pmf_k;
    R less  = R(1) - equal - more;
    if (less < R(0)) less = R(0);

    return { less, equal, more };
  }

}
