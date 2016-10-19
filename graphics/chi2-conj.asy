import graph;
import math;

real X=10;

unitsize(6*2.54cm/(X));
defaultpen(fontsize(12pt));

int k=4;
real x=5;
real y=0.5367762;
real z=0.1026063;

real f(real x) { return (x^(k/2-1)*exp(-x/2))/(2^(k/2)*gamma(k/2)); }

real scale=4/f(k-2);

real g(real x) { return scale*f(x); }
real c(real x) { return scale*z; }

path xax=(0,0)--(X,0);
draw(xax,Arrows);

yaxis("$y$",Arrows);

draw(graph(g,0,X,operator ..));

draw((y,0)--(y,scale*z),dashed);
draw((x,0)--(x,scale*z),dashed);
draw((y,scale*z)--(x,scale*z),dashed);

label("$x$",(x,0),S);
label("$x^*$",(y,0),S);
