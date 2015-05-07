function L=mnapproxluck(x,mu,Sigma)
  [dim,nsamps]=size(x);
  sigma=chol(Sigma)';
  z=sigma\(x-mu*ones(1,nsamps));
  R=sqrt(sum(z.^2,'r'));
  L=0.5*(1+erf(R-sqrt(dim-0.5)));
endfunction
