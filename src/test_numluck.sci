exec("mnluck.sci",2);
exec("mnsamp.sci",2);
exec("mnprobln.sci",2);
exec("numlucksetup.sci",2);
exec("numluck.sci",2);

function ok=test_mnluck()
  nsamps=1000;

  mu=[0.8;-0.5];
  Sigma=[4 3;3 5];

  x=mnsamp(nsamps,mu,Sigma);
  luck=mnluck(x,mu,Sigma);

  problns=mnprobln(x,mu,Sigma);
  setup=numlucksetup(problns);
  nluck=numluck(problns,setup);

  z=(nluck-luck) ./ sqrt(luck.*(1-luck)./nsamps);

  assert_checktrue(max(abs(z)) < 5);

  ok=%T;
endfunction
