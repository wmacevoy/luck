clear();

exec("mulsamp.sci",2);
exec("mulprobln.sci",2);
exec("mulluck.sci",2);
exec("mulapproxluck.sci",2);
exec("numlucksetup.sci",2);
exec("numluck.sci",2);

// get sample set
nsamps=4;
ntrials=100;
p=[0.4;0.3;0.2;0.1];
x=mulsamp(nsamps,ntrials,p);
problns=mulprobln(x,p);
// setup for luck estimates
setup=numlucksetup(problns);

// estimate luck
L=mulluck(x,p);
sd=sqrt(L .* (1-L) ./ nsamps);
L1=numluck(problns,setup);
z1=(L1-L) ./ sd;
L2=mulapproxluck(x,p);
z2=(L2-L) ./ sd;

L3=L1+max(-1,min(1,(L2-L1) ./ sd)) .* sd;
z3=(L3-L) ./ sd;
dz=(L2-L1) ./ sd;

c=linspace(-3,3,101);
clf();
histplot(c,z',style=1);