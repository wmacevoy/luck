d=100; // dimension
n=10000; // number of samples 

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

luck=hypermat([2,2,n]);
luck2=hypermat([2,2,n]);
for i=1:n
    r=zeros(2,2);
    r(1,1)=norm(a1*(Z1(:,i)-mu1));
    r(1,2)=norm(a2*(Z1(:,i)-mu2));
    r(2,1)=norm(a1*(Z2(:,i)-mu1));
    r(2,2)=norm(a2*(Z2(:,i)-mu2));
    L=cdfgam("PQ",r.^2/2.0,d/2.0*ones(2,2),ones(2,2));
    z=r-sqrt(d-0.5);
    el=0.5*(1+erf(z));
        
    luck(:,:,i)=L;
    luck2(:,:,i)=el;
//    luck11(i)=L(1,1);
//    luck12(i)=L(1,2);
//    luck21(i)=L(2,1);
//    luck22(i)=L(2,2);
//    printf("%12.10f & %12.10f & %12.10f & %12.10f \\\\\n",luck11(i),luck12(i),luck21(i),luck22(i));
end

luck=max(0,luck);
luck=min(1,luck);
luck2=max(0,luck2);
luck2=min(1,luck2);


xy=["x","y"];
clf();
c = linspace(0,1.0,101);
for i=1:2
    for j=1:2
        subplot(2,2,2*(i-1)+j);
        histplot(c,matrix(luck(i,j,:),n,1),style=1);
//        histplot(c,matrix(luck2(i,j,:),n,1),style=2);
        xtitle("L^"+xy(j)+"(" + xy(i) + ")");

    end
end
