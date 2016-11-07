clear();
exec("mnsamp.sci",-1);
exec("zluck.sci",-1);

df=5;
nsamps=100;
population=100;
generations=100;

x=zeros(df,population);
zl=zeros(1,nsamps);
Sigma0=zeros(df,df);
mu0=zeros(df,1);
Sigma1=zeros(df,df,nsamps);
mu1=zeros(df,nsamps);
Sigma2=zeros(df,df,nsamps);
mu2=zeros(df,nsamps);
q=zeros(nsamps+1);
mates=zeros(2,nsamps+1);

// goal PDF
mu0=rand(df,1);
Sigma0=rand(df,df);
Sigma0=Sigma0'*Sigma0;

// random initial sample
mu1=rand(df,nsamps);
Sigma1=rand(df,df,nsamps);
for i=1:nsamps
  Sigma1(:,:,i)=Sigma1(:,:,i)'*Sigma1(:,:,i);
end

mu1(:,3)=mu0;
Sigma1(:,:,3)=Sigma0;

for generation=1:generations
  // sample for generation
  for i=1:nsamps
    x=mnsamp(population,mu1(:,i),Sigma1(:,:,i));
    zl(i)=norm(zluck(x,mu0,Sigma0)+sqrt(df-1/2))-sqrt(population*df-1/2);
  end
  // overall luck in goal distribution
  zltotal=norm(zl+sqrt(population*df-1/2))-sqrt(nsamps*population*df-1/2);
  zlbar=(1/sqrt(nsamps))*(zltotal+sqrt(nsamps*population*df-1/2))-sqrt(population*df-1/2);

  // probability of mating
  q(2:nsamps+1)=0.5*erfc(abs(zl)-min(abs(zl)));
  q=q./sum(q);
  // mate cumulative probability map
  mates=[cumsum(q);0:nsamps];
  // mate pairs + random errors
  for i=1:nsamps
    p1=ceil(interpln(mates,rand()));
    p2=ceil(interpln(mates,rand()));
    w=zeros(1,3);
    w(1)=q(p1+1);
    w(2)=q(p2+1);
    z2=sqrt((zl(p1)+sqrt(df-1/2)).^2 + (zl(p2)+sqrt(df-1/2)).^2)-sqrt(2*df-1/2);
    w(3)=0.5*erfc(zlbar-z2);
    w=w./sum(w);
    tmp=rand(df,1);
    mu2(:,i)=w(1)*mu1(:,p1)+w(2)*mu1(:,p2)+w(3)*tmp;
    tmp=rand(df,df);
    tmp=tmp'*tmp;
    Sigma2(:,:,i)=w(1)*Sigma1(:,:,p1)+w(2)*Sigma1(:,:,p2)+w(3)*tmp;
  end // mate pairs

  // update generation
  mu1=mu2;
  Sigma1=Sigma2;
  disp(zlbar);
end // generations



