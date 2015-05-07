clear();

exec("mulsamp.sci",2);
exec("mulprobln.sci",2);
exec("mulluck.sci",2);
exec("numlucksetup.sci",2);
exec("numluck.sci",2);

// get sample set
nsamps=10000;
ntrials=100;
p=[0.1;0.2;0.3;0.4];
x=mulsamp(nsamps,ntrials,p);
problns=mulprobln(x,p);

// setup for luck estimates
setup=numlucksetup(problns);

// estimate luck numerically
nluck=numluck(problns,setup);
nsd=sqrt(nluck .* (1-nluck) ./ nsamps);

// optional exact luck
luck=mulluck(x,p);
sd=sqrt(luck .* (1-luck) ./ nsamps);

// z is approximately normally disributed
z=(nluck-luck) ./ sd;

c=linspace(-3,3,101);
clf();
histplot(c,z',style=1);