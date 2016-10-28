exec("chi2cdf.sci",-1);
exec("chi2conj.sci",-1);
exec("nonzero.sci",-1);

function [L,U]=chi2luck(x,k)
  y=chi2conj(x,k);
  [a,b]=chi2cdf(x,k);
  [A,B]=chi2cdf(y,k);
  L=abs(A-a);
  U=(A+b).*bool2s(x>=y)+(a+B).*bool2s(x<y);
endfunction
