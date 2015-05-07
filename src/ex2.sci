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

