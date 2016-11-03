exec("chi2luck.sci",-1);

k=4;
xmin=(max(0,sqrt(k-2)-3))^2;
xmax=(sqrt(k-2)+3)^2;

x=[xmin:0.1:xmax];

[L,U]=chi2luck(x,k);
[q,p]=chi2cdf(x,k);
P=exp(chi2probln(x,k));
clf();
plot(x,L,'k-',x,P,'k--',x,p,'k:');
xs2pdf(0,"../img/chi2.pdf");
