exec("zluck.sci",-1);

function [L,U]=normalluck(x,mu,Sigma)
  zl=zluck(x,mu,Sigma);
  L=0.5*erfc(-zl);
  U=0.5*erfc(zl);
endfunction
