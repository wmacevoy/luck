#pragma once

#include "luck.h"
#include <cmath>

namespace luck {

  // ---------------------------------------------------------------------------
  // Cauchy distribution with location x0 and scale gamma > 0:
  //   pdf(x) = 1 / (pi * gamma * (1 + ((x - x0)/gamma)^2))
  //   CDF(x) = 1/2 + (1/pi) * atan((x - x0)/gamma)
  // ---------------------------------------------------------------------------
  template <typename R>
  R cauchy_log_pdf(R x0, R gamma, R y) {
    R z = (y - x0) / gamma;
    return -std::log(R(M_PI) * gamma) - std::log(R(1) + z * z);
  }

  template <typename R>
  R cauchy_pdf(R x0, R gamma, R y) {
    R z = (y - x0) / gamma;
    return R(1) / (R(M_PI) * gamma * (R(1) + z * z));
  }

  // ---------------------------------------------------------------------------
  // Cauchy luck mass at observation y of Cauchy(x0, gamma):
  //   less  = P(pdf(X) <  pdf(y)) = 1 - (2/pi) * atan(|y - x0|/gamma)
  //   equal = 0                                  (continuous: |Omega| = 0)
  //   more  = P(pdf(X) >  pdf(y)) =     (2/pi) * atan(|y - x0|/gamma)
  //
  // The pdf is unimodal at x0 and strictly decreasing in |x - x0|, so the band
  // {x: pdf(x) >= pdf(y)} is the symmetric interval [x0 - d, x0 + d] with
  // d = |y - x0|. Its mass is CDF(x0+d) - CDF(x0-d) = (2/pi)*atan(d/gamma).
  // ---------------------------------------------------------------------------
  template <typename R>
  mass<R> cauchy_mass(R x0, R gamma, R y) {
    if (gamma <= R(0)) return { R(0), R(0), R(0) };

    R d = std::abs(y - x0) / gamma;
    R more = R(2) / R(M_PI) * std::atan(d);
    R less = R(1) - more;
    return { less, R(0), more };
  }

}
