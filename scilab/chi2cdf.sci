function [P,Q]=chi2cdf(x,k)
  [dim,nsamps]=size(x);
  one=ones(1,nsamps);
  [P,Q]=cdfgam("PQ",x./2,k/2*one,one);
endfunction
