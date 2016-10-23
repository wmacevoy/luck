clear();

exec("chi2luck.sci",-1);

k=4;
clf();
xmin=(max(0,sqrt(k-2)-3))^2;
xmax=(sqrt(k-2)+3)^2;

x=[xmin:0.01:xmax];

[L,U]=chi2luck(x,k);
[q,p]=chi2cdf(x,k);
P=exp(chi2probln(x,k));
plot(x',[L',P',p']);
