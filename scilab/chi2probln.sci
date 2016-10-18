exec("nonzero.sci",2);
function lnp=chi2probln(x,k)
  lnp=log(nonzero(x,1e-300)).*(k/2-1)-x./2-((k/2)*log(2)+gammaln(k/2));
endfunction
