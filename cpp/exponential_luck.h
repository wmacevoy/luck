#pragma once

#include "luck.h"
#include <cmath>

namespace luck {

  // ---------------------------------------------------------------------------
  // Exponential distribution -- the continuous analog of the geometric, with
  // rate parameter lambda > 0 on x >= 0:
  //   pdf(x) = lambda * exp(-lambda * x)
  // ---------------------------------------------------------------------------
  template <typename R>
  R exponential_log_pdf(R lambda, R y) {
    return std::log(lambda) - lambda * y;
  }

  template <typename R>
  R exponential_pdf(R lambda, R y) {
    return lambda * std::exp(-lambda * y);
  }

  // ---------------------------------------------------------------------------
  // Exponential luck mass at observation y of Exp(lambda):
  //   less  = P(pdf(X) <  pdf(y)) = P(X > y) = exp(-lambda*y)
  //   equal = P(pdf(X) == pdf(y)) = 0           (continuous: |Omega| = 0)
  //   more  = P(pdf(X) >  pdf(y)) = P(X < y) = 1 - exp(-lambda*y)
  //
  // The pdf is strictly decreasing on the support, so the "band" {x: pdf(x) >=
  // pdf(y)} is [0, y] and the pieces are just the CDF and survival evaluated at
  // y -- no walk, no twin, no equal mass.
  // ---------------------------------------------------------------------------
  template <typename R>
  mass<R> exponential_mass(R lambda, R y) {
    if (lambda <= R(0)) return { R(0), R(0), R(0) };
    if (y < R(0))       return { R(0), R(0), R(0) };  // outside support

    R surv = std::exp(-lambda * y);   // P(X > y)
    R cdf  = R(1) - surv;             // P(X < y)
    return { surv, R(0), cdf };
  }

}
