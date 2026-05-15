#include "exponential_luck.h"
#include <cmath>
#include <cstdio>

using luck::mass;
using luck::exponential_mass;
using luck::exponential_pdf;

// Brute-force reference: midpoint-rule integration of the pdf over a generous
// support window, bucketing each slice by pdf(x) vs pdf(y). For a continuous
// distribution `equal` is measure-zero and only picks up integration error.
template <typename R>
mass<R> ref_exponential_mass(R lambda, R y) {
  R pdfy = exponential_pdf<R>(lambda, y);
  R less = 0, equal = 0, more = 0;
  R xmax = R(60) / lambda;   // exp(-60) ~ 1e-26 tail
  int    N  = 2'000'000;
  R dx = xmax / R(N);
  R tol = R(1e-10);
  for (int i = 0; i < N; ++i) {
    R x = (R(i) + R(0.5)) * dx;
    R pdfx = exponential_pdf<R>(lambda, x);
    R prob = pdfx * dx;
    R d = pdfx - pdfy;
    if      (std::abs(d) <= tol * pdfy) equal += prob;
    else if (d < 0)                     less  += prob;
    else                                more  += prob;
  }
  return { less, equal, more };
}

static int fails = 0;
template <typename R>
void check(const char *label, R lambda, R y) {
  auto a = exponential_mass<R>(lambda, y);
  auto b = ref_exponential_mass<R>(lambda, y);
  // The closed form is exact; the reference is limited by the midpoint-rule
  // slice that straddles y, whose error is ~ pdf(y) * dx ~ 1.5e-5 here.
  R tol = R(1e-4);
  bool ok = std::abs(a.less  - b.less)  < tol
         && std::abs(a.equal - b.equal) < tol
         && std::abs(a.more  - b.more)  < tol;
  std::printf("%-22s lambda=%6.3f y=%7.3f  less=%.10f equal=%.10f more=%.10f  %s\n",
              label, (double)lambda, (double)y,
              (double)a.less, (double)a.equal, (double)a.more,
              ok ? "ok" : "FAIL");
  if (!ok) {
    std::printf("  ref:                                  less=%.10f equal=%.10f more=%.10f\n",
                (double)b.less, (double)b.equal, (double)b.more);
    ++fails;
  }
}

int main() {
  // lambda = 1: y at the median (ln 2) splits 50/50.
  check<double>("rate1-median",   1.0, std::log(2.0));
  check<double>("rate1-mode",     1.0, 0.0);
  check<double>("rate1-y0.5",     1.0, 0.5);
  check<double>("rate1-y2",       1.0, 2.0);
  check<double>("rate1-y5",       1.0, 5.0);

  // Other rates.
  check<double>("rate2-median",   2.0, std::log(2.0) / 2.0);
  check<double>("rate2-y1",       2.0, 1.0);
  check<double>("rate0.5-y1",     0.5, 1.0);
  check<double>("rate0.5-y10",    0.5, 10.0);
  check<double>("rate0.1-y5",     0.1, 5.0);
  check<double>("rate0.1-y30",    0.1, 30.0);

  if (fails == 0) {
    std::printf("\nall checks passed\n");
    return 0;
  } else {
    std::printf("\n%d FAILURES\n", fails);
    return 1;
  }
}
