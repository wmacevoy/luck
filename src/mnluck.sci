function L=mnluck(x,mu,Sigma)
  [dim,nsamps]=size(x);
  sigma=chol(Sigma)';
  z=sigma\(x-mu*ones(1,nsamps));
  R=sqrt(sum(z.^2,'r'));
  L=cdfgam("PQ",R.^2/2,dim/2*ones(1,nsamps),ones(1,nsamps));
endfunction
