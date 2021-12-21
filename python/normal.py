from sympy import *
from sympy.utilities.lambdify import *
from sympy.functions.special.gamma_functions import *
from sympy.functions.special.error_functions import *

r=Symbol('r',real=True)
df=Symbol('df',positive=True,integer=True)
luck_multivariate_normal = lambdify([r,df],lowergamma(df/2, r**2/2)/gamma(df/2))
luck_approximate_multivariate_normal = lambdify([r,df],(1/2)*(1+erf(r-sqrt(df-1/2))))

print("df,min,max,err,approx_err,ratio")
for ndf in [floor(exp(3)),floor(exp(4)),floor(exp(5))]:
    spread = 3
    samples = 201
    nrs = [sqrt(ndf-0.5)+spread*i/((samples-1)/2) for i in range(-((samples-1)//2),(samples-1)//2+1)]

    Lmvn = [N(luck_multivariate_normal(nr,ndf)) for nr in nrs]
    Lamvn = [N(luck_approximate_multivariate_normal(nr,ndf)) for nr in nrs]
    err = [abs(Lmvn[i]-Lamvn[i]) for i in range(len(Lmvn))]
    max_err = max(err)
    approx_err = N(1/(20*sqrt(ndf)))

    print(str(ndf) + "," + str(min(Lmvn)) + "," + str(max(Lmvn)) + "," + str(max_err) + "," + str(approx_err) + "," + str(max_err/approx_err))
