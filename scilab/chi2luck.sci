exec("chi2probln.sci",2);
exec("chi2cdf.sci",2);
exec("nonzero.sci",2);

function [L,U]=chi2luck(x,k)
  [dim,nsamps]=size(x);
  lnc=chi2probln(x,k);
  lny=(1/(k/2-1))*(chi2probln(x,k)+(k/2)*log(2)+gammaln(k/2));
  y=max(exp(lny),2*(k-2)-x);
  for i=1:10
    y = y + (2*y ./ nonzero(y-k+2,1e-6)) .* (chi2probln(y,k)-lnc);
  end
  [a,b]=chi2cdf(x,k);
  [A,B]=chi2cdf(y,k);
  L=abs(B-b);
  U=(A+b).*bool2s(x>=y)+(a+B).*bool2s(x<y);
endfunction
