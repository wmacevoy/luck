exec("chi2cdf.sci",-1);
exec("chi2luck.sci",-1);
exec("mulsamp.sci",-1);

pbad=[0.8,0.9,0.9,0.9,0.9,1.0]';
pbad=pbad ./ sum(pbad);

nsamps=3;
nrolls=45;

rolls=mulsamp(nsamps,nrolls,pbad);

pgood=[1,1,1,1,1,1]';
pgood=pgood./sum(pgood);

erolls=nrolls*(p*ones(1,nsamps));

x0=sum(((rolls-erolls) .^2) ./ erolls, 'r');

df=5;
[qval0,pval0]=chi2cdf(x0,df);
[L0,U0]=chi2luck(x0,df);

rolls1 = ((nrolls-pmodulo(nrolls,6))/6 + ([0:5] < pmodulo(nrolls,6)))';
rolls1 = rolls1 * ones(1,nsamps);
x1=sum(((rolls1-erolls) .^2) ./ erolls, 'r');
[qval1,pval1]=chi2cdf(x1,df);
[L1,U1]=chi2luck(x1,df);

