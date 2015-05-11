exec("mulluckrec.sci",2);

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
