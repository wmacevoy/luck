function luck=numluck(problns,setup)
  luck=max(0,min(1,interpln(setup,-problns)));
endfunction
