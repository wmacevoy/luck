clear();

exec("chi2luck.sci",-1);

k=4;
x=[0:0.01:10];
[L,U]=chi2luck(x,k);
[q,p]=chi2cdf(x,k);
P=exp(chi2probln(x,k));
plot(x',[L',P',p']);
