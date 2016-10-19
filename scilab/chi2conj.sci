exec("chi2probln.sci",-1);
exec("chi2cdf.sci",-1);
exec("nonzero.sci",-1);

function y=chi2conj(x,k)
  lnp=chi2probln(x,k);
  y=exp((1/(k/2-1))*(lnp+(k/2)*log(2)+gammaln(k/2)));
  y=max(y,2*(k-2)-x);
  iterate=%T;
  i=0;
  cutoff=sqrt(%eps);
  while (iterate)
    z=y;
    dln=2*z ./ nonzero(z-k+2,%eps);
    y = z +  dln .* (chi2probln(z,k)-lnp);
    i=i+1;
    iterate = (i<4) & (max(abs(y-z))>cutoff);
  end
endfunction
