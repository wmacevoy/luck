exec("zluck.sci",2);

function [L,U]=mnluck(x,mu,Sigma)
  [df,nsamps]=size(x);
  one=ones(1,nsamps);
  zl=zluck(x,mu,Sigma);
  r2=(zl+sqrt(df/2)).^2;
  [L,U]=cdfgam("PQ",r2/2,(df/2)*one,one);
endfunction
