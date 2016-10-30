mu=[1;1.1];
sigma=[1;1.1];
n=1000;
x=grand(n,'mn',mu(1),sigma(1)^2);
y=grand(n,'mn',mu(2),sigma(2)^2);

zl=zeros(2,2);

zl(1,1) = sqrt(sum(((x-mu(1)) ./ sigma(1)).^2))-sqrt(n-0.5);
zl(1,2) = sqrt(sum(((x-mu(2)) ./ sigma(2)).^2))-sqrt(n-0.5);
zl(2,1) = sqrt(sum(((y-mu(1)) ./ sigma(1)).^2))-sqrt(n-0.5);
zl(2,2) = sqrt(sum(((y-mu(2)) ./ sigma(2)).^2))-sqrt(n-0.5);

nl = 0.5*(1+erf(zl));



