#include "neg_binomial_luck.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using luck::mass;
using luck::neg_binomial_mass;
using luck::neg_binomial_pmf;

// Brute-force reference: enumerate far into the tail and bucket each term
// by comparison with pmf(k). jmax is taken comfortably past mean + 20*sd
// (and at least 300 to handle small-p slow decays).
template <typename R>
mass<R> ref_neg_binomial_mass(int r, R p, R q, int k) {
  R pk = neg_binomial_pmf<R>(r, p, q, k);
  R less = 0, equal = 0, more = 0;
  R tol = R(64) * std::numeric_limits<R>::epsilon();
  R mean = R(r) * q / p;
  R sd   = std::sqrt(R(r) * q) / p;
  int jmax = std::max(300, int(std::ceil(mean + R(20) * sd)) + 60);
  if (jmax < k + 4) jmax = k + 4;
  for (int j = 0; j <= jmax; ++j) {
    R pj = neg_binomial_pmf<R>(r, p, q, j);
    R d = pj - pk;
    if      (std::abs(d) <= tol * pk) equal += pj;
    else if (d < 0)                   less  += pj;
    else                              more  += pj;
  }
  return { less, equal, more };
}

static int fails = 0;
template <typename R>
void check(const char *label, int r, R p, int k) {
  R q = R(1) - p;
  auto a = neg_binomial_mass<R>(r, p, q, k);
  auto b = ref_neg_binomial_mass<R>(r, p, q, k);
  R tol = R(1e-11);
  bool ok = std::abs(a.less  - b.less)  < tol
         && std::abs(a.equal - b.equal) < tol
         && std::abs(a.more  - b.more)  < tol;
  std::printf("%-22s r=%3d p=%.4f k=%4d  less=%.10f equal=%.10f more=%.10f  %s\n",
              label, r, (double)p, k,
              (double)a.less, (double)a.equal, (double)a.more,
              ok ? "ok" : "FAIL");
  if (!ok) {
    std::printf("  ref:                                  less=%.10f equal=%.10f more=%.10f\n",
                (double)b.less, (double)b.equal, (double)b.more);
    ++fails;
  }
}

int main() {
  // r=1 reduces to geometric.
  check<double>("geom-k0",       1, 0.5, 0);
  check<double>("geom-k1",       1, 0.5, 1);
  check<double>("geom-k5",       1, 0.5, 5);

  // r=2, p=0.5: mode_real = 0 -> tied peaks at 0 and 1.
  check<double>("r2.5-peak-low",  2, 0.5, 0);
  check<double>("r2.5-peak-high", 2, 0.5, 1);
  check<double>("r2.5-k2",        2, 0.5, 2);
  check<double>("r2.5-k5",        2, 0.5, 5);
  check<double>("r2.5-k15",       2, 0.5, 15);

  // r=5, p=0.3: mode at 9 (non-integer mode_real ~ 8.33), no tie.
  check<double>("r5.3-mode",      5, 0.3, 9);
  check<double>("r5.3-below",     5, 0.3, 3);
  check<double>("r5.3-above",     5, 0.3, 20);
  check<double>("r5.3-k0",        5, 0.3, 0);

  // r=10, p=0.5: mode_real = 8 (integer) -> tied peaks at 8 and 9.
  check<double>("r10.5-peak-low",  10, 0.5, 8);
  check<double>("r10.5-peak-high", 10, 0.5, 9);
  check<double>("r10.5-below",     10, 0.5, 4);
  check<double>("r10.5-above",     10, 0.5, 18);

  // r=3, p=0.9: mode_real < 0, mode at 0.
  check<double>("r3.9-k0", 3, 0.9, 0);
  check<double>("r3.9-k1", 3, 0.9, 1);
  check<double>("r3.9-k4", 3, 0.9, 4);

  // r=20, p=0.4: larger distribution.
  check<double>("r20.4-mode",  20, 0.4, 28);
  check<double>("r20.4-below", 20, 0.4, 10);
  check<double>("r20.4-above", 20, 0.4, 60);

  if (fails == 0) {
    std::printf("\nall checks passed\n");
    return 0;
  } else {
    std::printf("\n%d FAILURES\n", fails);
    return 1;
  }
}
