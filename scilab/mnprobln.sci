function y=mnprobln(x,mu,Sigma)
  [dim,nsamps]=size(x);
  sigma=chol(Sigma)';
  z=sigma\(x-mu*ones(1,nsamps));
  R2=sum(z.^2,'r');
  y=-R2/2-(dim/2*log(2*%pi)+sum(log(diag(sigma))));
endfunction
