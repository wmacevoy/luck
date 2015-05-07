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
