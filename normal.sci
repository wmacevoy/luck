d=25; // dimension
n=1000; // number of samples 

// random mean & covariance
mu1=rand(d,1);
sigma1=rand(d,d);
sigma1=sigma1'*sigma1;
a1=chol(inv(sigma1));

mu2=rand(d,1);
sigma2=rand(d,d);
sigma2=sigma2'*sigma2;
a2=chol(inv(sigma2));

Z1=grand(n, "mn", mu1,sigma1);
Z2=grand(n, "mn", mu2,sigma2);

// luckij is the luck of seeing a sample generated using distribution
// i as an outcome in distribution j.
for i=1:n
    z11(i)=norm(a1*(Z1(:,i)-mu1))-sqrt(d-0.5);
    luck11(i)=0.5*(1+erf(z11(i)));
    z12(i)=norm(a2*(Z1(:,i)-mu2))-sqrt(d-0.5);
    luck12(i)=0.5*(1+erf(z12(i)));
    z21(i)=norm(a1*(Z2(:,i)-mu1))-sqrt(d-0.5);
    luck21(i)=0.5*(1+erf(z21(i)));
    z22(i)=norm(a2*(Z2(:,i)-mu2))-sqrt(d-0.5);
    luck22(i)=0.5*(1+erf(z22(i)));
//    printf("%12.10f & %12.10f & %12.10f & %12.10f \\\\\n",luck11(i),luck12(i),luck21(i),luck22(i));
end

c = linspace(-0.0001,1.0001,21);
f1=scf(1);
clf();
histplot(c,luck11',style=1);
xtitle("L^x(x)");
f2=scf(2);
clf();
histplot(c,luck12',style=1);
xtitle("L^x(y)");
f3=scf(3)
clf();
histplot(c,luck21',style=1);
xtitle("L^y(x)");
f4=scf(4);
clf();
histplot(c,luck22',style=1);
xtitle("L^y(y)");


