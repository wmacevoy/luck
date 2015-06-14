function ok=test_mul(verbose)
  if ~exists("verbose","local") then
    verbose=%F
  end
  test_mulprobln(verbose=verbose);
  test_mulsamp(verbose=verbose);
  test_mulluck();
  ok=%T;
endfunction
