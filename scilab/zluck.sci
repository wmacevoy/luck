function [zl,df,r]=zluck(x,mu,Sigma)
  [df,nsamps]=size(x);
  sigma=chol(Sigma)';
  one=ones(1,nsamps);
  z=sigma\(x-mu*one);
  r=sqrt(sum(z.^2,'r'));
  zl=r-sqrt(df-1/2);
endfunction
