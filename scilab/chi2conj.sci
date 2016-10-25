exec("chi2probln.sci",-1);
exec("chi2cdf.sci",-1);
exec("nonzero.sci",-1);

function y=chi2conj(x,k)
  [df,nsamps]=size(x);
  lnp=chi2probln(x,k);
  y0=exp((1/(k/2-1))*(lnp+(k/2)*log(2)+gammaln(k/2)));
  y0=max(y0,max(0,2*sqrt(k-2)-sqrt(x)) .^ 2);
  y=y0;
  iterate=%T;
  i=0;
  cutoff=sqrt(%eps);
  while (iterate)
    z=y;
    df=chi2probln(z,k)-lnp;
    dln=2*z ./ nonzero(z-k+2,cutoff);
    y = z +  dln .* df;
    i=i+1;
    iterate = (i<1000) & or(abs(df)>cutoff);
  end
endfunction
