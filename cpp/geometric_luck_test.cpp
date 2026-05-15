#include "geometric_luck.h"
#include <cmath>
#include <cstdio>
#include <limits>

using luck::mass;
using luck::geometric_mass;

// Brute-force reference: sum the pmf far enough that the remaining tail is
// below the comparison tolerance, bucketing each term against pmf(k).
template <typename R>
mass<R> ref_geometric_mass(R p, int k) {
  R q  = R(1) - p;
  R pk = p * std::pow(q, R(k));
  R less = 0, equal = 0, more = 0;
  R tol = R(64) * std::numeric_limits<R>::epsilon();
  // Tail q^j < 1e-20 when j > -20*ln(10)/ln(q); pad generously.
  int jmax = (q > R(0)) ? int(std::ceil(R(46) / -std::log(q))) + 60 : 1;
  if (jmax < k + 4) jmax = k + 4;
  for (int j = 0; j <= jmax; ++j) {
    R pj = p * std::pow(q, R(j));
    R d = pj - pk;
    if      (std::abs(d) <= tol * pk) equal += pj;
    else if (d < 0)                   less  += pj;
    else                              more  += pj;
  }
  return { less, equal, more };
}

static int fails = 0;
template <typename R>
void check(const char *label, R p, int k) {
  auto a = geometric_mass<R>(p, k);
  auto b = ref_geometric_mass<R>(p, k);
  R tol = R(1e-12);
  bool ok = std::abs(a.less  - b.less)  < tol
         && std::abs(a.equal - b.equal) < tol
         && std::abs(a.more  - b.more)  < tol;
  std::printf("%-18s p=%7.4f k=%4d  less=%.10f equal=%.10f more=%.10f  %s\n",
              label, (double)p, k,
              (double)a.less, (double)a.equal, (double)a.more,
              ok ? "ok" : "FAIL");
  if (!ok) {
    std::printf("  ref:                          less=%.10f equal=%.10f more=%.10f\n",
                (double)b.less, (double)b.equal, (double)b.more);
    ++fails;
  }
}

int main() {
  // Fair-ish coin
  check<double>("half-k0",   0.5, 0);
  check<double>("half-k1",   0.5, 1);
  check<double>("half-k5",   0.5, 5);
  check<double>("half-k20",  0.5, 20);

  // Asymmetric
  check<double>("p=.3-k0",   0.3, 0);
  check<double>("p=.3-k3",   0.3, 3);
  check<double>("p=.3-k10",  0.3, 10);

  // Rare success
  check<double>("p=.1-k0",   0.1, 0);
  check<double>("p=.1-k5",   0.1, 5);
  check<double>("p=.1-k50",  0.1, 50);
  check<double>("p=.01-k0",  0.01, 0);
  check<double>("p=.01-k99", 0.01, 99);

  // Likely success
  check<double>("p=.9-k0",   0.9, 0);
  check<double>("p=.9-k1",   0.9, 1);
  check<double>("p=.9-k5",   0.9, 5);

  // Edge: certain success
  check<double>("p=1-k0",    1.0, 0);
  check<double>("p=1-k3",    1.0, 3);

  if (fails == 0) {
    std::printf("\nall checks passed\n");
    return 0;
  } else {
    std::printf("\n%d FAILURES\n", fails);
    return 1;
  }
}
