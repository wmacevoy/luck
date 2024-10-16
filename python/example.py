import math

def factorial(n):
    ans = 1
    for i in range(2,n+1):
        ans = ans*i
    return ans

def combinations(n, k):
    return factorial(n)//(factorial(k)*factorial(n-k))

def binomial(n, k, p):
    return combinations(n,k)*p**k*(1.0-p)**(n-k)

X = [0,1,2,3,4,5,6,7,8]

def prob(x) :
    return binomial(8,x,0.5)

def prob_Omega(x):
    global X
    ans = 0.0
    px = prob(x)
    for y in X:
        py = prob(y)
        if px < py:
            ans += py
    return ans
            

def prob_omega(x):
    global X
    ans = 0.0
    ans = 0.0
    px = prob(x)
    for y in X:
        py = prob(y)
        if px == py:
            ans += py
    return ans

def luck(x):
    return prob_Omega(x) + 0.5*prob_omega(x)

def coins(x):
    L = luck(x)
    U = 1.0 - L
    return math.log2(L/U)


for x in X:
    print(f"x={x}, p(x)={prob(x)}, L(x)={luck(x)}, L_C(x)={coins(x)}")

