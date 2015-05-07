test_mul();
nsamps=10000;
ntrials=100;
p=[0.1;0.2;0.3;0.4];
verbose=%T;

x=mulsamp(nsamps,ntrials,p);
problns=mulprobln(x,p);
setup=numlucksetup(problns);
  
x0=[13;]
probllns0=mulprobln(x0,p);
nlucks=numluck(problns,setup);
nsd=sqrt(nlucks .* (1-nlucks) ./ nsamps);

lucks=mulluck(x,p);
sd=sqrt(lucks .* (1-lucks) ./ nsamps);
z=(nlucks-lucks) ./ sd;

  if verbose then
    printf("nluck=%0.15f, nsd=%0.15f, luck=%0.15f, sd=%0.15f, z=%15f\n",nlucks,nsd,lucks,sd,z);
  end
endfunction
