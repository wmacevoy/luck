exec("mnluck.sci",2);
exec("mnprobln.sci",2);

function ok=test_mnluck()
  mu=0.5;
  Sigma=3.3;

  x=[-4:4];

  L1=mnluck(x,mu,Sigma);
  L2=erf(abs((x-mu)/(2*Sigma)));

  assert_checkalmostequal(L1,L2);

  mu=[1;-1];
  Sigma=[5 6;6 9];

  n1=101;
  y1=linspace(-10,10,n1);
  n2=101;
  y2=linspace(-10,10,n2);
  y=zeros(2,n1*n2);
  for i1=1:n1
    b=(i1-1)*n2;
    y(1,b+1:b+n2)=y1(i1);
    y(2,b+1:b+n2)=y2;
  end
  dy=((y1(n1)-y1(1))/(n1-1))*((y2(n2)-y2(1))/(n2-1));
  p=exp(mnprobln(y,mu,Sigma))*dy;

  for k = [3000,4000,5000,6000]
    L1=mnluck(y(:,k),mu,Sigma);
    L2=sum(p.* bool2s(p > p(k)));
    assert_checkalmostequal(L1,L2,0.01);
  end

  ok=%T;
endfunction