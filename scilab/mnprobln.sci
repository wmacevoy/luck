function y=mnprobln(x,mu,Sigma)
  [df,nsamps]=size(x);
  sigma=chol(Sigma)';
  z=sigma\(x-mu*ones(1,nsamps));
  r2=sum(z.^2,'r');
  y=-r2/2-(df/2*log(2*%pi)+sum(log(diag(sigma))));
endfunction
