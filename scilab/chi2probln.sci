exec("nonzero.sci",-1);
function lnp=chi2probln(x,k)
  lnx=log(nonzero(x,number_properties('tiny')));
  lnp=lnx.*(k/2-1)-x./2-((k/2)*log(2)+gammaln(k/2));
endfunction
