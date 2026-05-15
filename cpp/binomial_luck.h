#pragma once

#include "luck.h"
#include <cmath>
#include <limits>

namespace luck {


  // ---------------------------------------------------------------------------
  // Continuous binomial log-pmf via lgamma. Defined for x in (-1, n+1).
  // ---------------------------------------------------------------------------
  template <typename R>
  R binomial_log_pmf(int n, R p, R q, R x) {
    return std::lgamma(R(n) + R(1))
         - std::lgamma(x + R(1))
         - std::lgamma(R(n) - x + R(1))
         + x         * std::log(p)
         + (R(n) - x) * std::log(q);
  }

  template <typename R>
  R binomial_pmf(int n, R p, R q, int k) {
    return std::exp(binomial_log_pmf<R>(n, p, q, R(k)));
  }

  // Sum of pmf(j) over the half-open index range [a, b).
  //
  // This is a self-contained O(b-a) fallback. The point of binomial_mass is the
  // *pattern*: less/equal/more are determined by a contiguous index band, so the
  // band sums are just CDF differences. A backend with a binomial CDF (= the
  // regularized incomplete beta, I_q(n-k,k+1)) should drop this enumeration and
  // evaluate the band endpoints directly -- O(1) in n instead of O(n).
  template <typename R>
  R binomial_sum(int n, R p, R q, int a, int b)
  {
    if (b-a < 4) {
      R s = R(0);
      for (int k=a; k<b; ++k) {
	s += binomial_pmf(n,p,q,k);
      }
      return s;
    } else {
      int c = (a+b)/2;
      return binomial_sum(n,p,q,a,c)+binomial_sum(n,p,q,c,b);
    }
  }

  // ---------------------------------------------------------------------------
  // Binomial luck mass at observation k of Bin(n,p):
  //   less  = sum_{j: pmf(j) <  pmf(k)} pmf(j)
  //   equal = sum_{j: pmf(j) == pmf(k)} pmf(j)   (= pmf(k) generically, 2*pmf(k) at a tied peak)
  //   more  = sum_{j: pmf(j) >  pmf(k)} pmf(j)
  // Strategy: the pmf is unimodal, so {j: pmf(j) >= pmf(k)} is a contiguous
  // band with k as one endpoint and its "twin" el (the equal-height index on
  // the far side of the mode) as the other. Find el on the (monotone) log-pmf
  // flank, then less/equal/more are fixed by that band.
  //
  // Here the band is summed with binomial_sum (O(n)). In a backend with a
  // binomial CDF this is the useful reference pattern: with the band endpoints
  // known, less = CDF(el-1) [+ tail], more = CDF(k-1) - CDF(el), etc. -- two
  // O(1) CDF calls, no enumeration. See the note on binomial_sum above.
  // ---------------------------------------------------------------------------
  template <typename R>
  mass<R> binomial_mass(int n, R p, R q,int k) {
    if (n < 0 || k < 0 || k > n) return { R(0), R(0), R(0) };
    if (p <= R(0)) return { (k == 0) ? R(0) : R(1), (k == 0) ? R(1) : R(0), R(0) };
    if (p >= R(1)) return { (k == n) ? R(0) : R(1), (k == n) ? R(1) : R(0), R(0) };

    int mode = std::floor(R(n+1)*p);
    bool mode_pair = ((mode > 0) && R(mode)*q == R(n-mode+1)*p);

    R log_pmf_k = binomial_log_pmf<R>(n, p, q, k);
    R pmf_k = std::exp(log_pmf_k);

    if (k == mode || (mode_pair && k == mode-1)) {
      R eq = (mode_pair) ? R(2)*pmf_k : pmf_k;
      return {R(1)-eq,eq,R(0)};
    }

    auto g = [&](R x) { return binomial_log_pmf<R>(n, p, q, x) - log_pmf_k; };
    // Tolerance on the log-pmf scale; scaled by n so it tracks the lgamma
    // cancellation error and still recognizes genuine integer twins.
    R eps = R(64)*R(n+1)*std::numeric_limits<R>::epsilon();

    // el is the farthest index on the opposite side of the mode from k that
    // still satisfies pmf(el) >= pmf(k). The high-or-equal band is [el,k] (or
    // [k,el]); g is monotone on each flank, so a real-valued warm start plus
    // single-step walks snap el onto the exact integer twin.
    int el;
    if (k > mode) {
      // el on the rising left flank: smallest j in [0,mode] with g(j) >= -eps.
      if (g(R(0)) >= -eps) {
	el = 0;
      } else {
	el = int(std::lround(interval_solve<R>(g, R(0), R(mode-int(mode_pair)), R(1.0/1024))));
	if (el < 0) el = 0; else if (el > mode) el = mode;
      }
      while (el-1 >= 0 && g(R(el-1)) >= -eps) { --el; }
      while (el < mode && g(R(el))   <  -eps) { ++el; }
    } else {
      // el on the falling right flank: largest j in [mode,n] with g(j) >= -eps.
      if (g(R(n)) >= -eps) {
	el = n;
      } else {
	el = int(std::lround(interval_solve<R>(g, R(mode), R(n), R(1.0/1024))));
	if (el < mode) el = mode; else if (el > n) el = n;
      }
      while (el+1 <= n && g(R(el+1)) >= -eps) { ++el; }
      while (el > mode && g(R(el))   <  -eps) { --el; }
    }

    // same: pmf(el) ties pmf(k), so el belongs to equal rather than more.
    int same = std::abs(g(R(el))) <= eps;

    mass<R> m{R(0),R(0),R(0)};
    if (el < k) {
      // less = [0,el-1] u [k+1,n];  more = [el(+1 if tie),k-1]
      m.less = binomial_sum(n,p,q,0,el) + binomial_sum(n,p,q,k+1,n+1);
      m.more = binomial_sum(n,p,q,el+same,k);
    } else {
      // less = [0,k-1] u [el+1,n];  more = [k+1,el(-1 if tie)]
      m.less = binomial_sum(n,p,q,0,k) + binomial_sum(n,p,q,el+1,n+1);
      m.more = binomial_sum(n,p,q,k+1,el+1-same);
    }
    m.equal = R(1+same)*pmf_k;

    return m;
  }

}
