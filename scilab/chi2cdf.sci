function [p,q]=chi2cdf(x,k)
  [dim,nsamps]=size(x);
  [p,q]=cdfgam("PQ",x./2,k/2*ones(1,nsamps),ones(1,nsamps));
endfunction
