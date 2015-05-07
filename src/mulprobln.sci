function y=mulprobln(x,p)
  [nprob,nsamp]=size(x);
  pln=log(p);
  y=zeros(1,nsamp);
  for i=1:nsamp
    xi=x(:,i);
    ni=sum(xi);
    y(i)=gammaln(ni+1)+sum(xi.*pln-gammaln(xi+1));
  end
endfunction
