exec("mnprobln.sci",2);

function ok=test_mnprobln()
  mu=[1;-1];
  Sigma=[5 6;6 9];

  n1=51;
  y1=linspace(-10,10,n1);
  n2=51;
  y2=linspace(-10,10,n2);
  y=zeros(2,n1*n2);
  for i1=1:n1
    b=(i1-1)*n2;
    y(1,b+1:b+n2)=y1(i1);
    y(2,b+1:b+n2)=y2;
  end
  dy=((y1(n1)-y1(1))/(n1-1))*((y2(n2)-y2(1))/(n2-1));

  p=exp(mnprobln(y,mu,Sigma))*dy;

  m0=sum(p);
  m1=zeros(2,1);
  m1(1)=sum(p .* y(1,:))/m0;
  m1(2)=sum(p .* y(2,:))/m0;

  m2=zeros(2,2);  
  for i=1:2
    for j=1:2
      m2(i,j)=sum(p .* (y(i,:)-mu(i)).*(y(j,:)-mu(j)))/m0;
    end
  end

  assert_checkalmostequal(m0,1,0.025);
  assert_checkalmostequal(m1,mu,0.025);
  assert_checkalmostequal(m2,Sigma,0.025);

  ok=%T;
endfunction