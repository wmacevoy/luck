#pragma once

#include "luck.h"
#include <cmath>
#include <limits>

namespace luck {

  // ---------------------------------------------------------------------------
  // Geometric distribution: number of failures k = 0, 1, 2, ... before the
  // first success in a sequence of Bernoulli(p) trials. pmf(k) = p * q^k with
  // q = 1 - p. Continuous form via x*log(q) + log(p), defined for x > -1.
  // ---------------------------------------------------------------------------
  template <typename R>
  R geometric_log_pmf(R p, R q, R x) {
    return x * std::log(q) + std::log(p);
  }

  template <typename R>
  R geometric_pmf(R p, R q, int k) {
    return p * std::pow(q, R(k));
  }

  // ---------------------------------------------------------------------------
  // Geometric luck mass at observation k of Geom(p):
  //   less  = sum_{j: pmf(j) <  pmf(k)} pmf(j) = q^{k+1}
  //   equal = sum_{j: pmf(j) == pmf(k)} pmf(j) = p q^k
  //   more  = sum_{j: pmf(j) >  pmf(k)} pmf(j) = 1 - q^k
  //
  // The pmf is strictly decreasing, so the "band" {j: pmf(j) >= pmf(k)} is just
  // [0,k], and each piece is closed form: more is the geometric partial sum
  // 1 - q^k, less is the geometric tail q^{k+1}. No walk, no twin search.
  // ---------------------------------------------------------------------------
  template <typename R>
  mass<R> geometric_mass(R p, int k) {
    if (k < 0)         return { R(0), R(0), R(0) };
    if (p <= R(0))     return { R(0), R(0), R(0) };           // no success ever
    if (p >= R(1)) {                                          // pmf(0)=1, else 0
      if (k == 0)      return { R(0), R(1), R(0) };
      return { R(0), R(0), R(1) };
    }

    R q     = R(1) - p;
    R qk    = std::pow(q, R(k));
    R equal = p * qk;
    R more  = R(1) - qk;     // sum_{j=0}^{k-1} p q^j
    R less  = q * qk;        // q^{k+1}

    return { less, equal, more };
  }

}
