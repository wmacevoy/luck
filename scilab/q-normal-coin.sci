p=[1/4, 3/4];
x=[0,1];
mu=p(1)*x(1)+p(2)*x(2);
sigma=sqrt(p(1)*(x(1)-mu)^2+p(2)*(x(2)-mu)^2);
printf("mu=%f sigma=%f\n",mu,sigma);
