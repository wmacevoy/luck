function y=mulprobln(x,p)
  [nprob,nsamp]=size(x);
  pln=log(p);
  y=zeros(1,nsamp);
  for i=1:nsamp
    xi=x(:,i);
    ni=sum(xi);
    y(i)=gammaln(ni+1)+sum(xi.*pln-gammaln(xi+1));
  end
endfunction

function ok=test_mulprobln(n,verbose)
  if ~exists("n","local") then
    n=4
  end
  if ~exists("verbose","local") then
    verbose=%F
  end

  ok=%F;
  p=[0.3;0.7];
  for x1=0:n
    for x2=0:n
      x=[x1;x2];
      c=gamma(x1+x2+1)/(gamma(x1+1)*gamma(x2+1));
      ans1=mulprobln(x,p);
      ans2=log(c*p(1)^x1*p(2)^x2);
      if verbose
        printf("x=[%d,%d] ans1=%f, ans2=%f\n",x(1),x(2),ans1,ans2);
      end
      assert_checkalmostequal(ans1,ans2,0,sqrt(%eps));
    end
  end

  p=[0.1;0.5;0.4];
  for x1=0:n
    for x2=0:n
      for x3=0:n
        x=[x1;x2;x3];
        c=gamma(x1+x2+x3+1)/(gamma(x1+1)*gamma(x2+1)*gamma(x3+1));
	ans1=mulprobln(x,p);
        ans2=log(c*p(1)^x1*p(2)^x2*p(3)^x3);
        if verbose
	  printf("x=[%d,%d,%d] ans1=%f, ans2=%f\n",x(1),x(2),x(3),ans1,ans2);
        end
        assert_checkalmostequal(ans1,ans2,0,sqrt(%eps));
      end
    end
  end
  ok=%T;
endfunction

function x=mulsamp(nsamps,ntrials,p)
  x=grand(nsamps,"mul",ntrials,p(1:(length(p)-1)));
endfunction

function test_mulsamp(verbose)
  p=[0.1;0.5;0.4];
  nsamp=1000000;
  ntrials=10000;
  x=mulsamp(nsamp,ntrials,p);  
  mean1=sum(x,'c') ./ nsamp;
  mean2=ntrials*p;

  assert_checkalmostequal(mean1,mean2,max(sqrt(1/nsamp),sqrt(1/ntrials)));
  
  cov1=zeros(length(p),length(p));
  cov2=zeros(length(p),length(p));

  for i=1:length(p)
    for j=1:length(p)
      cov1(i,j)= sum((x(i,:)-mean2(i)).*(x(j,:)-mean2(j)))/nsamp;
      if i == j then
        cov2(i,j)=ntrials*p(i)*(1-p(i));
      else
        cov2(i,j)=-ntrials*p(i)*p(j);
      end
    end
  end

  assert_checkalmostequal(cov1,cov2,max(sqrt(1/nsamp),sqrt(1/ntrials)));
endfunction

function el=mulluckrec(nb,mb,lpx,k,p,y,eps)
  [n,m]=size(lpx);
  el=zeros(n,m);
  for yk=0:nb-mb
     y(k)=yk;
     if k < length(p)-1 then
       el=el+mulluckrec(nb,mb+yk,lpx,k+1,p,y,eps);
     else
       y(length(p))=nb-mb-yk;
       lpy=mulprobln(y,p);
       c=0.5*bool2s(lpy > lpx-eps)+0.5*bool2s(lpy > lpx+eps);
       el = el + c .* exp(lpy);
     end
  end
endfunction

function el=mulluck(x,p,eps)
  if ~exists("eps","local") then
    eps=sqrt(%eps);
  end
  [nprobs,nsamps]=size(x);
  ntrials=sum(x,'r');
  min_ntrials=min(ntrials);
  max_ntrials=max(ntrials);
  assert_checkequal(min_ntrials,max_ntrials);
  el=mulluckrec(min_ntrials,0,mulprobln(x,p),1,p,zeros(length(p),1),eps);
endfunction

function ok=test_mulluck()
  ok=%F;
  x=[1;5;7;5];
  p=[0.1;0.2;0.3;0.4];
  ans1=mulluck(x,p);
  ans2=0.396540127837419;
  assert_checkalmostequal(ans1,ans2);

  x=[4 3
     4 5];

  p=[0.5
     0.5];

  ans1=mulluck(x,p);
  ans2=[(0.5)*(0.5^8*(gamma(8+1)/(gamma(4+1)*gamma(4+1)))),...
        (0.5^8)*(gamma(8+1)/(gamma(4+1)*gamma(4+1)))+...
        (0.5^8)*(gamma(8+1)/(gamma(5+1)*gamma(3+1)))];
  assert_checkalmostequal(ans1,ans2);

  ok=%T;
endfunction

function setup=numlucksetup(problns,eps)
  if ~exists("eps","local") then
    eps=sqrt(%eps);
  end

  n=length(problns);
  problns=gsort(problns);

  i=1;
  count=0;
  while i<=n
    j=i;
    while (j <= n-1 & abs(problns(j+1)-problns(i)) < eps)  
      j=j+1;
    end
    i=j+1;
    count=count+1;
  end

  setup=zeros(3,count);

  i=1;
  count=0;
  while i<=n
    j=i;
    while (j <= n-1 & abs(problns(j+1)-problns(i)) < eps)  
      j=j+1;
    end

    count=count+1;
    setup(1,count)=0.5*(problns(i)+problns(j));
    setup(2,count)=(i-1)/n;
    setup(3,count)=(j-i+1)/n;
    i=j+1;
  end
endfunction

function [lucks,abs_Omegas,abs_omegas]=numluck(problns,setup,eps)
  if ~exists("eps","local") then
    eps=sqrt(%eps);
  end

  [n,m]=size(problns);
  [three,nsetup]=size(setup);

  abs_Omega=zeros(n,m);
  abs_omega=zeros(n,m);

  for i=1:n
    for j=1:m
      lo=1;
      hi=nsetup;
      while %T
        mid=floor((hi+lo)/2)
        if setup(1,mid) > problns(i,j)+eps then
   	  if mid == lo then
            break
          end
          lo=mid;
        else
          if mid == hi then
            break
          end
          hi=mid;
        end
      end

      if abs(problns(i,j)-setup(1,hi))<eps then
        abs_Omega(i,j)=setup(2,hi);
        abs_omega(i,j)=setup(3,hi);
      elseif abs(problns(i,j)-setup(1,lo))<eps then
        abs_Omega(i,j)=setup(2,lo);
        abs_omega(i,j)=setup(3,lo);
      else
        abs_Omega(i,j)=setup(2,lo);
        abs_omega(i,j)=0.5*exp(problns(i,j));
      end
    end
  end
  lucks=abs_Omega+0.5*abs_omega;
endfunction

function ok=multest(verbose)
  if ~exists("verbose","local") then
    verbose=%F
  end
  test_mulprobln(verbose=verbose);
  test_mulsamp(verbose=verbose);
  test_mulluck();
  ok=%T;
endfunction

function [x,luck,z]=mulmain(nsamps,ntrials,p, verbose)
  multest();
  if ~exists("nsamps","local") then
    nsamps=10000;
  end

  if ~exists("ntrials","local") then
    ntrials=100;
  end

  if ~exists("p","local") then
    p=[0.1;0.2;0.3;0.4];
  end

  if ~exists("verbose","local") then
    verbose=%F;
  end

  x=mulsamp(nsamps,ntrials,p);
  problns=mulprobln(x,p);
  setup=numlucksetup(problns);
  
  x=[13;]
  nlucks=numluck(problns,setup);
  nsd=sqrt(nlucks .* (1-nlucks) ./ nsamps);

  lucks=mulluck(x,p);
  sd=sqrt(lucks .* (1-lucks) ./ nsamps);
  z=(nlucks-lucks) ./ sd;

  if verbose then
    printf("nluck=%0.15f, nsd=%0.15f, luck=%0.15f, sd=%0.15f, z=%15f\n",nlucks,nsd,lucks,sd,z);
  end
endfunction

mulmain(verbose=%T);

