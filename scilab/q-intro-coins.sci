// possible x's (domain), three trials with possible max's of 0 or 1

outcomes=[[0,0,0];[1,0,0];[0,1,0];[0,0,1];[0,1,1];[1,0,1];[1,1,0];[1,1,1]];

[cases,trials]=size(outcomes);

// compute the probability of each outcome
probs = ones(cases,1);

probPerTrial = [1/4; 3/4];

for trial = 1:trials
    probs = probs .* probPerTrial(outcomes(:,trial)+1);
end

// sort probabilities and outcomes in decreasing order of probability
[probs,ord]=gsort(probs,'r');
outcomes=outcomes(ord,:);

// compute luck, with a small numerical tolerance for equality in probability
Omega = zeros(cases,1);
omega = zeros(cases,1);
luck = zeros(cases,1);
relerr = 1e-9;

total = 0;
i = 1;
while i <= cases
    part = probs(i);
    pMin = (1-relerr)*probs(i);
    j = i + 1;
    while j <= cases && probs(j) >= pMin
        part = part + probs(j);
        j = j + 1;
    end
    Omega(i:j-1)=total;
    omega(i:j-1)=part;
    luck(i:j-1)=Omega(i)+0.5*omega(i);
    total = total + part;
    i = j;
end

printf("x,p(x),|Omega(x)|,|omega(x)|,L(x)\n");
for i = 1:cases
    printf("%d %d %d,%0.4f,%0.4f,%0.4f,%0.4f\n",outcomes(i,1),outcomes(i,2),outcomes(i,3),probs(i),Omega(i),omega(i),luck(i));
end
