// ns=[10,100,1000,10000,100000];
log10n = linspace(1,5,1000);
ns = floor(10 .^ log10n);
as = zeros(1,length(ns));
bs = zeros(1,length(ns));


for i = [1:length(ns)]
    n=ns(i);
    R0=sqrt(n-1/2);
    z=linspace(max(-R0,-10),10,1000);
    R=z+R0;
    L1=cdfgam("PQ",(R .^ 2) ./ 2, ones(1,length(R)) .* (n/2), ones(1,length(R)));
    mu=sqrt(n-1/2);
    f1=erf(R-mu);
    f2=1+f1;
    f3=(f2) ./ 2;
    L2=f3;
    mu=sqrt(n);
    f1=erf(R-mu);
    f2=1+f1;
    f3=(f2) ./ 2;
    L3=f3;
    
    as(i)=max(abs(L1-L2));
    bs(i)=max(abs(L1-L3));
end
fig=scf(0);
clf();
plot2d(log10n,[as ./ bs]);
xtitle("abs err ratio vs log_10 n")
xs2png(gcf(),'/home/user/projects/luck/img/estluck.png');

