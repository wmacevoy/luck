
exec("chi2luck.sci",-1);

function chi2plots(k)
clf();
xmin=(max(0,sqrt(k-2)-3))^2;
xmax=(sqrt(k-2)+3)^2;

x=[xmin:0.01:xmax];
[L,U]=chi2luck(x,k);
[q,p]=chi2cdf(x,k);
P=exp(chi2probln(x,k));
Lest = erf(abs(sqrt(x)-sqrt(k-2)));
plot(x',[L',Lest']);
disp(max(abs(L-Lest)));
endfunction
