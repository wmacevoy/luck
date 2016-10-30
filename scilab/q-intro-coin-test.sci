
// 'random' numbers
c=[0 1 1 0 0 1];

L=[127/128, 117/128, 81/128, 27/128];
a=c(1:2:length(c)-1);
b=c(2:2:length(c));
m=max(a,b);
s=sum(m);
printf("luck(%d,%d,%d)=%f\n",m(1),m(2),m(3),L(s+1));

