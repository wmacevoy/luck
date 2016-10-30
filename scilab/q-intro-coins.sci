p1=zeros(2,1);
p1(1)=1/4; // given
p1(2)=3/4; // given
p3=zeros(2,2,2);
for out1=1:2
  for out2=1:2
    for out3=1:2
      p3(out1,out2,out3)=p1(out1)*p1(out2)*p1(out3);
    end
  end
end

for out1=1:2
  for out2=1:2
    for out3=1:2
      L=0;
      for other1=1:2
        for other2=1:2
          for other3=1:2
            if (p3(other1,other2,other3) > p3(out1,out2,out3)) then
               L = L + p3(other1,other2,other3);
            elseif (p3(other1,other2,other3) == p3(out1,out2,out3)) then
	       L = L + 0.5*p3(other1,other2,other3);
            end
          end
        end
      end
      printf("L(%d,%d,%d)=%f\n",out1-1,out2-1,out3-1,L);
    end
  end
end
