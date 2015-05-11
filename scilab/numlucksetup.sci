function setup=numlucksetup(problns,eps)
  if ~exists("eps","local") then
    eps=sqrt(%eps);
  end

  n=length(problns);
  problns=gsort(problns);

  for pass=1:2
    if pass == 2 then
      setup=zeros(2,count)
    end
    i=1;
    count=0;
    while i<=n
      j=i;
      while (j <= n-1 & abs(problns(j+1)-problns(i)) < eps)  
        j=j+1;
      end
      count=count+1;
      if pass == 2 then
        setup(1,count)= -0.5*(problns(i)+problns(j));
        setup(2,count)=(i-1)/n+0.5*(j-i+1)/n;
      end
      i=j+1;
    end
  end
endfunction
