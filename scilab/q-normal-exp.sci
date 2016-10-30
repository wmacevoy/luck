p=[1/4, 3/4];
x=[0,1];
mu=p(1)*x(1)+p(2)*x(2);
sigma=sqrt(p(1)*(x(1)-mu)^2+p(2)*(x(2)-mu)^2);

a=[1 0 1 1 0 1 1 0 1 0 0 0 1 1 0 1 1 0 0 0 1 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1 0 0 1 ]';
b=[1 0 1 0 1 0 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 0 0 1 0 0 1 0 1 1 0 1 0 0 1 0 0 ]';
n=length(a);
m=max(a,b);
zl=(m-mu) ./ sigma - sqrt(0.5);
for i=1:n
  printf("zl(%d)=%f\n",i,zl(i));
end
zl_sum=sqrt(sum((zl+sqrt(0.5)) .^ 2))-sqrt(n-0.5);
printf("zl(total)=%f\n",zl_sum);
L_sum = 0.5*(1+erf(zl_sum));
printf("total normal luck=%f\n",L_sum);

