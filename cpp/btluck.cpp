#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <thread>

#define DEBUG 0

struct ln_pq {
  double ln_p,ln_q;
};
std::ostream& operator<<(std::ostream &out, const ln_pq &ln_pq) {
  return out << "{\"ln_p\":" << ln_pq.ln_p << ",\"ln_q\":" << ln_pq.ln_q << "}";
}

struct mass {
  double less,equal,more;
};

std::ostream& operator<<(std::ostream &out, const mass &m) {
  return out << "{\"less\":" << m.less << ",\"equal\":" << m.equal << ",\"more\":" << m.more << "}";
}

const double LN_0 = log(0.0);
const double LN_1 = log(1.0);

mass walk(double ln_p0,double ln_p1,double epsilon,
	   const std::vector<ln_pq>::const_iterator &a,
	   const std::vector<ln_pq>::const_iterator &b) {
  double ln_p_max = ln_p1;
  double ln_p_min = ln_p1;  
  for (auto i=a; i != b; ++i) {
    ln_p_max += i->ln_q;
    ln_p_min += i->ln_p;
  }
  if (DEBUG) {
    std::cout << ln_p_min << "<=ln_p0=" << ln_p0 << "<=" << ln_p_max << std::endl;
  }
  if (ln_p0 > ln_p_max+epsilon) {
    return {exp(ln_p1),0,0};
  }
  if (ln_p0 < ln_p_min-epsilon) {
    return {0,0,exp(ln_p1)};
  }
  if (ln_p_max-ln_p_min <= 2*epsilon) {
    return {0,exp(ln_p1),0};
  }
  
  mass m0 = walk(ln_p0,a->ln_p+ln_p1,epsilon,a+1,b);
  mass m1 = walk(ln_p0,a->ln_q+ln_p1,epsilon,a+1,b);
  
  return {m0.less + m1.less,m0.equal+m1.equal,m0.more+m1.more};
}

template<typename T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &v) {
    for (int i=0; i<v.size(); ++i) {
      out << ((i == 0) ? "[" : ",") << v[i];
    }
    return out << "]";
}

// p(k)=1/(k+y0) probability of max after y years,
// starting from y0 in the record
// input y0 100010101101 -- sequence of record years starting with y9

int main(int argc, const char *argv[])
{
  double y0 = atoi(argv[1]);
  for (int i = 2; i<argc; ++i) {
    std::string record = argv[i];
    double ln_p0 = 0;
    std::vector<ln_pq> ln_probs;

    for (int i=record.size()-1; i>=0; --i) {
      double y = i+y0;
      double p = 1/y;
      std::cout << "y=" << y << ",p=" << p << std::endl;
      if (y == 0) continue;
      double ln_p = log(p);
      double ln_q = log1p(-p);
      ln_p0 += (record[i] == '1') ? ln_p : ln_q;
      ln_probs.push_back({ln_p,ln_q});
    }
    std::cout << "p0 = " << exp(ln_p0) << std::endl;
    std::cout << "ln_probs" << ln_probs << std::endl;
    mass m = walk(ln_p0,0.0,pow(2.0,-30),ln_probs.begin(),ln_probs.end());

    std::cout << "mass(" << record << ")=" << m << std::endl;

  }
  return 0;
}
