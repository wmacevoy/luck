#include "cauchy_luck.h"
#include <cmath>
#include <cstdio>

using luck::mass;
using luck::cauchy_mass;

// Brute-force reference: integrate via the t = atan((x-x0)/gamma) change of
// variable so the Cauchy measure is uniform on (-pi/2, pi/2). Each slice
// contributes 1/N of the mass exactly -- no tail truncation, only the slice
// that straddles the boundary +/- d contributes O(1/N) error.
template <typename R>
mass<R> ref_cauchy_mass(R x0, R gamma, R y) {
  R d = std::abs(y - x0) / gamma;
  R less = 0, equal = 0, more = 0;
  int N = 2000000;
  R pi = std::acos(R(-1));
  R prob = R(1) / R(N);
  for (int i = 0; i < N; ++i) {
    R t = (R(i) + R(0.5)) / R(N) * pi - pi / R(2);  // t in (-pi/2, pi/2)
    R xd = std::abs(std::tan(t));                   // |x - x0|/gamma
    if      (xd < d) more += prob;                  // closer to mode -> higher pdf
    else if (xd > d) less += prob;
    else             equal += prob;                 // measure zero in practice
  }
  return { less, equal, more };
}

static int fails = 0;
template <typename R>
void check(const char *label, R x0, R gamma, R y) {
  auto a = cauchy_mass<R>(x0, gamma, y);
  auto b = ref_cauchy_mass<R>(x0, gamma, y);
  // Closed form is exact; reference is limited by the straddle slice ~ 1/N.
  R tol = R(1e-4);
  bool ok = std::abs(a.less  - b.less)  < tol
         && std::abs(a.equal - b.equal) < tol
         && std::abs(a.more  - b.more)  < tol;
  std::printf("%-22s x0=%6.2f g=%6.3f y=%8.3f  less=%.10f equal=%.10f more=%.10f  %s\n",
              label, (double)x0, (double)gamma, (double)y,
              (double)a.less, (double)a.equal, (double)a.more,
              ok ? "ok" : "FAIL");
  if (!ok) {
    std::printf("  ref:                                          less=%.10f equal=%.10f more=%.10f\n",
                (double)b.less, (double)b.equal, (double)b.more);
    ++fails;
  }
}

int main() {
  // Standard Cauchy: x0 = 0, gamma = 1. y = 1 is the quartile boundary
  // (atan(1) = pi/4, so more = 0.5).
  check<double>("std-mode",      0.0, 1.0,  0.0);
  check<double>("std-quartile",  0.0, 1.0,  1.0);
  check<double>("std-quartile-",0.0, 1.0, -1.0);
  check<double>("std-y3",        0.0, 1.0,  3.0);
  check<double>("std-y10",       0.0, 1.0, 10.0);
  check<double>("std-y100",      0.0, 1.0, 100.0);

  // Shifted / scaled.
  check<double>("shift-mode",    5.0, 2.0, 5.0);
  check<double>("shift-q",       5.0, 2.0, 7.0);    // y - x0 = gamma
  check<double>("shift-q-",      5.0, 2.0, 3.0);
  check<double>("shift-far",     5.0, 2.0, 25.0);

  // Narrow / wide scale.
  check<double>("narrow-y1",     0.0, 0.1,  1.0);
  check<double>("wide-y1",       0.0, 10.0, 1.0);

  if (fails == 0) {
    std::printf("\nall checks passed\n");
    return 0;
  } else {
    std::printf("\n%d FAILURES\n", fails);
    return 1;
  }
}
