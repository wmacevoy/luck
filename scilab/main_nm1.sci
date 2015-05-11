clear();

exec("mnluck.sci",2);
exec("mnapproxluck.sci",2);

for dim=[1,2,5,10,20,50]
  mu=zeros(dim,1);
  Sigma=eye(dim,dim);

  nsamps=10000;
  x=zeros(dim,nsamps);
  x(1,:)=linspace(0,floor(2*sqrt(dim)),nsamps);

  L1=mnluck(x,mu,Sigma);
  L2=mnapproxluck(x,mu,Sigma);

  plot(x(1,:)',[L1' L2']);
end
