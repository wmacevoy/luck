import graph;
import math;

real xmin=0,xmax=3;
real ymin=0,ymax=1;

unitsize(6*2.54cm/(xmax-xmin));
defaultpen(fontsize(12pt));

typedef real r_r(real x);
typedef bool3 r_b3(real x);

void labpath(pair lp, real x0, r_r g, string tex, pair dir)
{
  real m=(g(x0+1e-6)-g(x0-1e-6))/(2e-6);
  pair u=(-m,1)/sqrt(m*m+1);
  pair p0=(x0,g(x0));
  pair dr=lp-p0;
  if (dr.x*u.x+dr.y*u.y < 0) {
    u=-u;
  }

  draw(lp..(lp-0.1*dir)..((x0,g(x0))+0.1*u)..(x0,g(x0)),arrow=Arrow(TeXHead));
  label(tex,lp,dir);
}

real f(real x) { return erf(abs(x/sqrt(2))); };
real g(real x) { return (1/2)*(1+erf(x-sqrt(0.5))); }

draw(graph(f,xmin,xmax,operator ..));
labpath((0.75,0.25),0.25,f,"${\tt erf} |x/\sqrt{2} |$ -- exact",E);

draw(graph(g,xmin,xmax,operator ..),dashed);
labpath((0.75,0.8),0.25,g,"$(1+{\tt erf}(x-\sqrt{1/2}))/2$ -- estimate",N);

xaxis("$x$",Bottom,LeftTicks(Label(fontsize(8pt)),new real[]{0,1,2,3}));
yaxis("$L$",Left,RightTicks(Label(fontsize(8pt)),new real[]{0,0.2,0.4,0.6,0.8,1.0}));
