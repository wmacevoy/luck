#include "poisson_luck.h"
#include <cmath>
#include <cstdio>
#include <limits>

using luck::mass;
using luck::poisson_mass;
using luck::poisson_pmf;

// Brute-force reference: enumerate the pmf far into the tail and bucket each
// term by comparison with pmf(k). Uses the same lgamma-based evaluator as
// poisson_mass so the comparison is consistent.
template <typename R>
mass<R> ref_poisson_mass(R lambda, int k) {
  R pk = poisson_pmf<R>(lambda, k);
  R less = 0, equal = 0, more = 0;
  R tol = R(64) * std::numeric_limits<R>::epsilon();
  int jmax = int(std::ceil(lambda + R(15) * std::sqrt(lambda + R(1)))) + 30;
  for (int j = 0; j <= jmax; ++j) {
    R pj = poisson_pmf<R>(lambda, j);
    R d = pj - pk;
    if      (std::abs(d) <= tol * pk) equal += pj;
    else if (d < 0)                   less  += pj;
    else                              more  += pj;
  }
  return { less, equal, more };
}

static int fails = 0;
template <typename R>
void check(const char *label, R lambda, int k) {
  auto a = poisson_mass<R>(lambda, k);
  auto b = ref_poisson_mass<R>(lambda, k);
  R tol = R(1e-12);
  bool ok = std::abs(a.less  - b.less)  < tol
         && std::abs(a.equal - b.equal) < tol
         && std::abs(a.more  - b.more)  < tol;
  std::printf("%-20s lambda=%8.3f k=%4d  less=%.10f equal=%.10f more=%.10f  %s\n",
              label, (double)lambda, k,
              (double)a.less, (double)a.equal, (double)a.more,
              ok ? "ok" : "FAIL");
  if (!ok) {
    std::printf("  ref:                                less=%.10f equal=%.10f more=%.10f\n",
                (double)b.less, (double)b.equal, (double)b.more);
    ++fails;
  }
}

int main() {
  // Integer lambda: tied modes at lambda-1 and lambda.
  check<double>("int-peak-low",  4.0, 3);
  check<double>("int-peak-high", 4.0, 4);
  check<double>("int-below",     4.0, 2);
  check<double>("int-above",     4.0, 6);
  check<double>("int-tail-0",    4.0, 0);
  check<double>("int-tail-far",  4.0, 12);

  // Non-integer lambda: single mode.
  check<double>("frac-mode",     4.3, 4);
  check<double>("frac-below",    4.3, 1);
  check<double>("frac-above",    4.3, 7);
  check<double>("frac-tail-0",   4.3, 0);
  check<double>("frac-tail-far", 4.3, 14);

  // Small lambda: mode at 0.
  check<double>("small-mode",    0.5, 0);
  check<double>("small-k1",      0.5, 1);
  check<double>("small-k3",      0.5, 3);

  // lambda = 1: tied modes 0 and 1.
  check<double>("one-k0", 1.0, 0);
  check<double>("one-k1", 1.0, 1);
  check<double>("one-k2", 1.0, 2);

  // Larger lambda.
  check<double>("big-mode",  100.0, 100);
  check<double>("big-below", 100.0,  80);
  check<double>("big-above", 100.0, 120);
  check<double>("big-tail",  100.0, 150);
  check<double>("mid-int",    25.0,  25);
  check<double>("mid-below",  25.0,  15);

  if (fails == 0) {
    std::printf("\nall checks passed\n");
    return 0;
  } else {
    std::printf("\n%d FAILURES\n", fails);
    return 1;
  }
}
