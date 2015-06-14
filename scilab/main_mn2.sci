exec("mnprobln.sci",2);
exec("mnluck.sci",2);
exec("mnsamp.sci",2);
exec("mnr.sci",2);

dim=25; // dimension
nsamps=1; // number of samples 

mu=hypermat([dim,2]);
Sigma=hypermat([dim,dim,2]);
x=hypermat([dim,nsamps,2]);

for i=1:2
  mu(:,i)=rand(dim,1);
  tmp=rand(dim,dim);
  Sigma(:,:,i)=tmp'*tmp;
  x(:,:,i)=mnsamp(nsamps,mu(:,i),Sigma(:,:,i));
end

z=zeros(2,2);

for i=1:2
  for j=1:2
    R=mnr(x(:,:,i),mu(:,j),Sigma(:,:,j));
    z(i,j)=norm(R)-sqrt(nsamps*dim-0.5);
  end
end

L=0.5*(1+erf(z));


