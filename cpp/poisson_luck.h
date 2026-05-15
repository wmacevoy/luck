#pragma once

#include "luck.h"
#include <cmath>
#include <limits>

namespace luck {

  // ---------------------------------------------------------------------------
  // Continuous Poisson log-pmf via lgamma. Defined for x > -1.
  // ---------------------------------------------------------------------------
  template <typename R>
  R poisson_log_pmf(R lambda, R x) {
    return -lambda + x * std::log(lambda) - std::lgamma(x + R(1));
  }

  template <typename R>
  R poisson_pmf(R lambda, int k) {
    return std::exp(poisson_log_pmf<R>(lambda, R(k)));
  }

  // ---------------------------------------------------------------------------
  // Poisson luck mass at observation k of Poisson(lambda):
  //   less  = sum_{j: pmf(j) <  pmf(k)} pmf(j)
  //   equal = sum_{j: pmf(j) == pmf(k)} pmf(j)   (= pmf(k), or 2*pmf(k) at an
  //                                               integer-lambda tied peak)
  //   more  = sum_{j: pmf(j) >  pmf(k)} pmf(j)
  //
  // The pmf is unimodal, so {j: pmf(j) >= pmf(k)} is a contiguous band with k as
  // one endpoint and its "twin" el (the equal-height index on the far side of
  // the mode) as the other. Unlike the binomial case there is no root solve:
  // the pmf ratio pmf(j+1)/pmf(j) = lambda/(j+1) lets us just walk the band from
  // k across the mode, accumulating `more` term by term, and stop when the pmf
  // falls back below pmf(k). `less` is then the complement 1 - equal - more, so
  // the infinite tail is never summed.
  // ---------------------------------------------------------------------------
  template <typename R>
  mass<R> poisson_mass(R lambda, int k) {
    if (k < 0) return { R(0), R(0), R(0) };
    if (lambda <= R(0)) {
      if (k == 0) return { R(0), R(1), R(0) };
      return { R(0), R(0), R(1) };
    }

    int  mode      = int(std::floor(lambda));
    bool mode_pair = (mode > 0 && R(mode) == lambda);  // lambda integer: modes mode-1, mode

    R pmf_k = poisson_pmf<R>(lambda, k);

    if (k == mode || (mode_pair && k == mode - 1)) {
      R eq = (mode_pair ? R(2) : R(1)) * pmf_k;
      return { R(1) - eq, eq, R(0) };
    }

    R eps    = R(64) * std::numeric_limits<R>::epsilon();
    R cutoff = pmf_k - eps * pmf_k;
    R more   = R(0);
    R pmf_el = R(0);   // pmf at the far band endpoint el

    if (k < mode) {
      // rising flank: walk up through the mode, pmf(j+1) = pmf(j)*lambda/(j+1).
      R pj = pmf_k;
      for (int j = k; ; ++j) {
        pj *= lambda / R(j + 1);          // pj = pmf(j+1)
        if (pj < cutoff) break;           // j+1 left the band; previous index is el
        pmf_el = pj;
        more  += pj;
      }
    } else {
      // falling flank: walk down through the mode, pmf(j-1) = pmf(j)*j/lambda.
      R pj = pmf_k;
      for (int j = k; j > 0; --j) {
        pj *= R(j) / lambda;              // pj = pmf(j-1)
        if (pj < cutoff) break;
        pmf_el = pj;
        more  += pj;
      }
    }

    // `more` currently runs through el itself; el belongs to `equal` if it ties
    // pmf(k), otherwise it is genuinely a `more` term and stays put.
    int same = std::abs(pmf_el - pmf_k) <= eps * pmf_k;
    if (same) more -= pmf_el;

    R equal = R(1 + same) * pmf_k;
    R less  = R(1) - equal - more;
    if (less < R(0)) less = R(0);        // guard tiny rounding overshoot

    return { less, equal, more };
  }

}
