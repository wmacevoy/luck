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
