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
