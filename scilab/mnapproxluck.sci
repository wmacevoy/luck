exec("mnr.sci",2);

function L=mnapproxluck(x,mu,Sigma)
  [dim,nsamps]=size(x);
  R=mnr(x,mu,Sigma);
  L=0.5*(1+erf(R-sqrt(dim-0.5)));
endfunction
