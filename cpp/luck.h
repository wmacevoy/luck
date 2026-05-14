#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace luck {

  template <typename R>
  struct mass {
    R less, equal, more;
  };

  template <typename R>
  struct LUC {
    R luck, unluck, coin;
  };

  template <typename R>
  mass<R> operator+(const mass<R> &a, const mass<R> &b) {
    return {a.less + b.less, a.equal + b.equal, a.more + b.more};
  }

  template <typename R>
  mass<R> operator*(const R &a, const mass<R> &b) {
    return {a * b.less, a * b.equal, a * b.more};
  }

  template <typename R>
  std::ostream& operator<<(std::ostream &out, const mass<R> &m) {
    return out << "{\"less\":"  << m.less
               << ",\"equal\":" << m.equal
               << ",\"more\":"  << m.more << "}";
  }

  template <typename R>
  R luck_of(const mass<R> &m) {
    return m.more + R(0.5) * m.equal;
  }

  template <typename R>
  R unluck_of(const mass<R> &m) {
    return m.less + R(0.5) * m.equal;
  }

  template <typename R>
  R coin_of(const mass<R> &m) {
    R heq = R(0.5) * m.equal;
    return std::log2(m.more + heq) - std::log2(m.less + heq);
  }

  template <typename R>
  LUC<R> luc_of(const mass<R> &m) {
    R heq    = R(0.5) * m.equal;
    R lk     = m.more + heq;
    R unlk   = m.less + heq;
    return { lk, unlk, std::log2(lk) - std::log2(unlk) };
  }

  // bisection / quadratic interval solver for finding single root in [a,b] with f(a)*f(b) <= 0
  template <typename R, typename F>
  R interval_solve(F &&f, R a, R b, R eps = R(0)) {
    if (eps <= R(0)) {
      eps = R(4) * (b-a)*std::numeric_limits<R>::epsilon();
    }

    R fa = f(a);
    R fb = f(b);

    for (;;) {
      R h = b-a;
      if (std::abs(h) <= eps) break;

      R c = R(0.5)*(a+b);
      R fc=f(c);

      // if monotonic, try quadratic estimator, otherwise bisect
      if (fa <= fc && fc <= fb || fa >= fc && fc >= fb) {

	// centered difference 1st and 2nd derivatives
	R fp = (fb - fa) / h;
	R fpp = R(4)*(fa + fb - R(2) * fc) / (h * h);

	// 2nd order root estimate (delta from c)
	R dh = -fc/fp;
	dh += -R(0.5)*fpp/fp*(dh*dh);

	c=c+dh;
	if (c < a || b < c) {
	  // revert to linear (guaranteed to be in interval)
	  c=R(0.5)*(a+b)-fc/fp;
	}

	fc=f(c);
	// check if within tolerance
	if (std::abs(fc)<=std::abs(fp*eps)) {
	  return c;
	}

	// estimate bound on root
	dh = -R(2)*fc/fp;
	if (dh < R(0)) { // [c+dh,c] should bound root
	  if (a < c+dh) {
	    b = a; fb = fa; // save for reset
	    a = c + dh;
	    fa = f(a);
	    if (fa*fc > R(0)) {
	      // no sign change - reset a
	      a = b;
	      fa=fb;
	    }
	  }
	  b = c;
	  fb = fc;
	} else { // [c,c+dh] should bound root
	  if (c+dh < b) {
	    a = b; fa = fb; // save for reset
	    b = c + dh;
	    fb = f(b);
	    if (fb*fc > R(0)) {
	      // no sign change - reset b
	      b = a;
	      fb = fa;
	    }
	  }
	  a = c;
	  fa = fc;
	}
      } else {
	// bisect
	if (fa*fc <= R(0)) {
	  b = c;
	  fb = fc;
	} else {
	  a = c;
	  fa = fc;
	}
      }
    }

    return R(0.5)*(a+b);
  }
}
