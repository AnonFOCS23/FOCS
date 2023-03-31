"""Code providing the following output:
******* BISECTION METHOD USED *******

       converged: True
           flag: 'converged'
 function_calls: 44
     iterations: 42
           root: 0.21123119692353295

**** ratio =  0.7234860332935461 ****

case = NOT MONOTONIC

* (a, b, p) =  (0.789, 1.24, 0.421) *

(lambda_*, mu_*, nu_*) =  (0.41518612252479703, 0.0022528733893470987, 0.21123119692353295)
"""
import numpy as np
from scipy.optimize import bisect #bisection method for root finding
from math import exp
from math import expm1 #exp(x)-1
from math import log
from sys import exit

def Emax(): # E max V
    return 1-b*expm1(-p)+a*exp(-p)

def l():#lambda_*
    return 1+log((1+(b-a)*p)/(1+b*p))/p

def m():#mu_*
    return 1-log(1+b*p)/p
    
def q(z): #q(x,y,z)=q(lambda,mu,nu), x,y fixed global
    return .5*y**2-.5*z**2+z+(1/p-y)*(1/p+b)*exp(p*(y-1))\
            +((1/p+b)*(z-x)-a/p)*exp(p*(z-1))\
            -(1/p+b-a)*exp(p*(z-x))/p+1/p+b

def dq(z): #q'(x,z)=q'(lambda,nu), x fixed global
    return 1-z+((1+b*p)*(z-x)-a)*exp(p*(z-1))

def condI():#LHS of Condition I
    return ((1+b*p)/(1+(b-a)*p))*log((1+b*p)/(1+(b-a)*p))

def condII():#LHS of Condition II
    return (2-p)*(b-a)

def condIII():#LHS of Condition III
    return (1-(2-p)*(b-a))/(1+p*(b-a))

def condIV():#LHS of Condition IV
    return 2+p*b*(1-p*(b-a))

def maximise_bisect(): #m(a,b,p)
    nu, r = bisect(dq, y, x, xtol = 1e-13, rtol = 1e-14, full_output = True)
    #nu = nu_* , r = convergence details
    print('\n******* BISECTION METHOD USED *******\n\n', r)
    return q(nu), (x,y,nu)

a = .789
b = 1.24
p = .421
x = l() #x = lambda_*
y = m() #y = mu_*
#check of all conditions
if condI()>a*p or condII()>=1 or condIII()>=x or condIV()<0 or log(1+p*b)>=p: 
    exit('Instance not feasible')

emax = Emax()

if dq(y) <= 0:
    maxq, pos = q(y), (x, y, y) #monotonic case
    case = 'monotonic'
else:
    maxq, pos = maximise_bisect() #nonmonotonic case
    case = 'nonmonotonic'
    
M = maxq/emax #M(a,b,p)

print('\n***** ratio = ',M,'******\n\n* (a, b, p) = ',(a, b, p),'*\n')
if case == 'monotonic':
    print('(lambda_*, mu_*, mu_*) = ', pos,'\n')
    print('case = MONOTONIC')
else:
    print('(lambda_*, mu_*, nu_*) = ', pos,'\n')
    print('case = NOT MONOTONIC')
    print('\n***** BISECTION ROOT APPROXIMATION ERROR ****\n')
    print('|root-true_root|<=xtol+rtol*|root|<=1e-13+.3*1e-14<1e-12')