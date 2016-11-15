function zlbar=zluckbar(x,mu,Sigma)
  [df,nsamps]=size(x);
  zl=zluck(x,mu,Sigma);
  zltotal=norm(zl+sqrt(df-1/2))-sqrt(nsamps*df-1/2);
  zlbar=zltotal/sqrt(df*nsamps);
endfunction
