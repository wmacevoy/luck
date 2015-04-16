import graph;
import math;

real xmin=0;
real xmax=1;

real ymin= -1;
real ymax= 1;

unitsize(6*2.54cm/(xmax-xmin));
defaultpen(fontsize(10pt));

int n=8;
real p=0.5;
real q=1-p;

real log_gamma(real x)
{
  return log(gamma(x));
}

real log_multinomial(int [] counts, real [] probs) {
  int total=0;
  real log_ans = 0;

  for (int k=0; k<counts.length; ++k) { total += counts[k]; }
  
  log_ans += log_gamma(total+1);

  for (int k=0; k<counts.length; ++k) {
    log_ans -= log_gamma(counts[k]+1);
    log_ans += counts[k]*log(probs[k]);
  }
  return log_ans;
}

real multinomial(int [] counts, real [] probs) {
  return exp(log_multinomial(counts,probs));
}

real eps=1e-6;

bool eq(real a, real b)
{
  return fabs(b-a)<eps;
}

real [] probs = {};
real [] cprobs = {};
real [] lucks = {};
int [] xs = {};
int [] ys;

for (int x=0; x<=n; ++x) {
  int [] multinomial_counts = {x,n-x};
  real [] multinomial_probs = {p,1-p};
  real p=multinomial(multinomial_counts,multinomial_probs);
  xs.push(x);
  probs.push(p);
  lucks.push(0);
  cprobs.push(0);
}

string fmt(real x,int n)
{
  if (x < 0) return "-" + fmt(-x,n);
  for (int i=0; i<n; ++i) x *= 10;
  string ans=string(floor(x+0.5));
  while (length(ans) < n+1) {
    ans = "0" + ans;
  }
  ans = insert(ans,length(ans)-n,".");
  return ans;
}

bool order(int x,int y) { return probs[x] > probs[y]; }

ys=sort(xs,order);

int i=0;
while (i < ys.length) {
  int j=i;
  while (j+1 < ys.length && eq(probs[ys[j+1]],probs[ys[i]])) { ++j; }
  real luck=(i>0?cprobs[ys[i-1]]:0)+probs[ys[i]]*(j-i+1)/2.0;
  while (i<=j) {
    cprobs[ys[i]]=(i>0?cprobs[ys[i-1]]:0)+probs[ys[i]];
    lucks[ys[i]]=luck;
    ++i;
  }
}

pair narrow(pair range,real h)
{
  return (min(range.y,range.x+h),max(range.x,range.y-h));
}

void interval(int k)
{
  pair range=(cprobs[ys[k]]-probs[ys[k]],cprobs[ys[k]]);
  real h=0.005;
  pair narrowed=narrow(range,h);
  real luck=lucks[ys[k]];
  real lucky=0.05;

  draw((narrowed.x,lucky*2)--(narrowed.y,lucky*2));
  if (range.y-range.x > 0.1) {
    label("$p="+fmt(probs[ys[k]],2)+"$",((range.x+range.y)/2,lucky*2),N);
  }
  if (k == 0 || luck != lucks[ys[k-1]]) {
    int w=0;
    while (k+w < ys.length && lucks[ys[k+w]] == lucks[ys[k]]) ++w;
    pair luck_range=(range.x,range.x+w*probs[ys[k]]);
    pair luck_narrowed=narrow(luck_range,h);
    draw((luck_narrowed.x,lucky*2-h)--(luck_narrowed.x,lucky)--(luck,lucky)--(luck_narrowed.y,lucky)--(luck_narrowed.y,lucky*2-h));

    draw((luck,lucky-h)--(luck,lucky+h));
    //    draw((0,-(lucks.length-k)*h)--(luck,-(lucks.length-k)*h));

    if (luck_range.y - luck_range.x > 0.1) {
      string L="L";
      L=L+"(x=";
      for (int el=0; el<w; ++el) {
	if (el > 0) L=L+"\mbox{ or }";
	L=L+string(ys[k+el]);
      }
      L=L+")";
      label("$"+L+"="+fmt(luck,2)+"$",(luck,0),N);
    } else if (luck_range.y - luck_range.x > 0.01) {
      label("$\cdots$",(luck,0),N);      
    }

    //    dot((luck,lucky));
    //    draw((luck,lucky)--(narrowed.x,lucky));
    //    draw((narrowed.x,lucky)--(narrowed.x,lucky-h));
    //    draw((luck,lucky)--(narrowed.y,lucky));
    //    draw((narrowed.y,lucky)--(narrowed.y,lucky-h));
  }
}

//draw((0,0)--(-0.25,0),arrow=Arrow(TeXHead));
//draw((1,0)--(1.25,0),arrow=Arrow(TeXHead));
for (int k=0; k<ys.length; ++k) {
  interval(k);
}


// dot((x+h,g2(x+h)));
