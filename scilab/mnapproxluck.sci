exec("zluck.sci",2);

function [L,U]=mnapproxluck(x,mu,Sigma)
  zl=zluck(x,mu,Sigma);
  L=0.5*erfc(-zl);
  U=0.5*erfc(zl);
endfunction
