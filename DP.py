"""Code providing the following output
*** Hard instance ***
n =  1000000 
a =  0.789 
b =  1.24 
p =  0.421

*** Acceptance times ***
k_1000000 =  2253
bar_k_1000000 =  211231
j_1000000 =  415187

***Gambler-to-prophet ratio***
0.723494757359636
"""
import numpy as np
from matplotlib import pyplot as plt

def Emax(a,b,p,n):
    return n*(1-np.power(1-1/n**2,n))+b*(np.power(1-1/n**2,n)\
            -np.power(1-p/n-1/n**2,n))+a*np.power(1-p/n-1/n**2,n)

def Eopt(a,b,p,n):#E gamma_1(X_1)
    return max(a,phi[0])/(n+1)+(1-1/(n+1))*(1/n+p*max(b,phi_[0])/n\
            +(1-p/n-1/n**2)*phi_[0]) 

def dp(a,b,p,n):#dynamic program
   phi=np.empty(n)
   phi_=np.empty(n)
   phi[n-1]=(1+b*p)/n
   phi_[n-1]=a
   for i in range(n-2,-1,-1):
       phi[i]=1/n+p*max(b,phi[i+1])/n+(1-p/n-1/n**2)*phi[i+1]
       phi_[i]=max(a,phi[i+1])/(n-i+1)+(1-1/(n-i+1))\
                *(1/n+p*max(b,phi_[i+1])/n+(1-p/n-1/n**2)*phi_[i+1])
   return phi,phi_

n=1000000
a=.789
b=1.24
p=.421

phi,phi_=dp(a,b,p,n) 
emax=Emax(a,b,p,n)
eopt=Eopt(a,b,p,n)

s=0
for i in phi:#find j_n
    if a>=i:
        j=s
        break
    s+=1
s=0   
for i in phi:#find k_n
    if  b>=i:
        k=s
        break
    s+=1
s=0
for i in phi_:#find bar_k_n
    if  b>=i:
        k_=s
        break
    s+=1
    
A=[a for i in range(n)]
B=[b for i in range(n)] 
plt.figure(dpi=1200)  #1200 for high res
plt.plot(phi_, color='orange', linestyle=':', linewidth=1) 
plt.plot(phi, color='blue', linestyle=':', linewidth=1)
plt.plot(B, color='red', linewidth=0.3, label='b')
plt.plot(A, color='green', linewidth=0.3, label='a')
plt.xlabel('k')
plt.legend(loc='lower left')

print('\n*** Hard instance ***')
print('n = ', n, '\na = ', a, '\nb = ', b, '\np = ', p)
print('\n*** Acceptance times ***')
print(f'k_{n} = ', k)
print(f'bar_k_{n} = ', k_)
print(f'j_{n} = ', j)
print('\n***Gambler-to-prophet ratio***')
print(eopt/emax)