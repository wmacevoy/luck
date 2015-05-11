function L=mulapproxluck(x,p)
  [nprobs,nsamps]=size(x);
  ntrials=sum(x,'r');
  mu=p*ntrials;
  z=(x-mu) ./ sqrt(mu);
  R2=sum(z.^2,'r');
  one=ones(1,nsamps);
  L=cdfgam("PQ",R2/2,((nprobs-1)/2)*one,one);
endfunction
