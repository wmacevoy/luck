exec("zluck.sci",2);

function [L,U]=mnluck(x,mu,Sigma)
  [df,nsamps]=size(x);
  one=ones(1,nsamps);
  [zl,df,r]=zluck(x,mu,Sigma);
  [L,U]=cdfgam("PQ",r.^2/2,(df/2)*one,one);
endfunction
