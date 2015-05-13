exec("mnr.sci",2);

function [L,U]=mnluck(x,mu,Sigma)
  [dim,nsamps]=size(x);
  R=mnr(x,mu,Sigma);
  [L,U]=cdfgam("PQ",R.^2/2,dim/2*ones(1,nsamps),ones(1,nsamps));
endfunction
