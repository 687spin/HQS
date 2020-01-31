# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:08:58 2019

@author: iwh
"""

from scipy.special import lpmn, clpmn
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import cmath
import time
from sympy import *
from mpmath import *

#Input constants
a = 1.0 #sphere radius
#H_0 = 1200.0 #Applied DC field, Oersted
M = 141.7 #Saturation magnetization, A/m
gamma = 2.8e6 #Hz / Gauss


n = 3
m = 0
s =n-m
ho = 2.2e4
sign = -1

def zero_finder(x,y):
    zeros = []
    for i in range(len(y) -1):
        if y[i] * y[i+1] < 0:
            zeros.append((x[i] + x[i+1])/2)
    return zeros

def y(F,n,m,H_0,sign):
    OmH = H_0/(4*np.pi*M) - (1/3)
    Omega = F + OmH
    kappa = -OmH/(F*(F+2*OmH))
    nu = -(F + OmH)/(F*(F+2*OmH))
    xi = np.sqrt(1 + 1/kappa)
    #print(1 + 1/kappa)
    nu = -(F + OmH)/(F*(F+2*OmH))
    return (n+1+sign*m*nu)*lpmn(m,n,xi)[0][m][n]+xi*lpmn(m,n,xi)[1][m][n]
    #return (n+1+sign*m*nu)*legenp(m,n,xi)+ xi * diff(assoc_legendre(n,m,xi),xi)
    #lpmn returns 2 arrays, lpmn[0] is the associated Legendre polynomoal, and lpmn[1] is the derivatives.
    #each array has all m,n combinations possible so its necessarry to specify lpmn[0][m][n] to get the desired value for the function
    
def y_t(F,n,m,H_0,sign):
    OmH = H_0/(4*np.pi*M) - (1/3)
    kappa = -OmH/(F*(F+2*OmH))
    nu = -(F + OmH)/(F*(F+2*OmH))
    xi = np.sqrt(abs(1 + 1/kappa))
    nu = -(F + OmH)/(F*(F+2*OmH))
    if kappa <= -1:
        return (n+1+sign*m*nu)*lpmn(m,n,xi)[0][m][n]+xi*lpmn(m,n,xi)[1][m][n]
    else:
        return -((n+1+sign*m*nu)*lpmn(m,n,xi)[0][m][n]+xi*lpmn(m,n,xi)[1][m][n])
    return 9 + (3*O)
    #lpmn returns 2 arrays, lpmn[0] is the associated Legendre polynomoal, and lpmn[1] is the derivatives.
    #each array has all m,n combinations possible so its necessarry to specify lpmn[0][m][n] to get the desired value for the function

def y_2(f,n,m,H_0,sign):
    Omega = f/(4*np.pi*gamma*M)
    #print ('Omega: {}'.format(Omega))
    H_i = H_0 - (4*np.pi/3)*M
    Omega_H = H_i/(4*np.pi*M)
    #print ('OmegaH: {}'.format(Omega_H))
    nu = Omega/(Omega_H**2 - Omega**2)
    kappa = Omega_H/(Omega_H**2 - Omega**2)
    xi=np.sqrt(1.0 + 1.0/kappa)
    #print (Omega - Omega_H)
    #return n + 1 + xi*lpmn(m,n,xi)[1][m][n]/lpmn(m,n,xi)[0][m][n] + sign*m*nu
    return (n+1+sign*m*nu)*lpmn(m,n,xi)[0][m][n]+xi*lpmn(m,n,xi)[1][m][n]
    #return (n+1+sign*m*nu)*lpmn(m,n,xi)[0][m][n]+xi*lpmn(m,n,xi)[1][m][n]
    #lpmn returns 2 arrays, lpmn[0] is the associated Legendre polynomoal, and lpmn[1] is the derivatives.
    #each array has all m,n combinations possible so its necessarry to specify lpmn[0][m][n] to get the desired value for the function
    

def y_3(f,n,m,H_0,sign):
    Omega = f*1000000000/(4*np.pi*gamma*M)
    H_i = H_0 - (4*np.pi/3)*M
    Omega_H = H_i/(4*np.pi*M)
    nu = Omega/(Omega_H**2-Omega**2)
    kappa = Omega_H/(Omega_H**2-Omega**2)
    xi=np.sqrt(1.0 + 1.0/kappa)
    f = lambda xi: legenp(n,m, xi)
    #return n + 1 + xi*lpmn(m,n,xi)[1][m][n]/lpmn(m,n,xi)[0][m][n] +sign*m*nu
    ft = lambda xi: (n+1+sign*m*nu)*f(n,m,xi) + xi*diff(f(n,m,xi), xi)
    print(f(xi))
    return f(xi)

    #lpmn returns 2 arrays, lpmn[0] is the associated Legendre polynomoal, and lpmn[1] is the derivatives.
    #each array has all m,n combinations possible so its necessarry to specify lpmn[0][m][n] to get the desired value for the function

def MSM_freq(H0, n , m, sign):
    if n - m < 2:
        Freqs = np.empty((len(H0), 1))
    else:
        Freqs = np.empty((len(H0), (n-m)//2 + 1))
    y_val = np.empty(len(x))
    for f, ho in enumerate(H0):
        for i, j in enumerate(x):
            y_val[i] = y(j, n, m, ho, sign)
        res = zero_finder(x, y_val)
        #print (res)
        #plt.plot(x, y_val)
        #plt.savefig(r'X:\Shared Data\YBCO\MSM_Freqencies\H0_{}_n_{}_m_{}.png'.format(round(ho,2), n, m), format = 'png')
        #plt.close()
        Freqs[f, :] = res
    return Freqs

#H0 = np.linspace(0.2*np.pi*M, 6*4*np.pi*M, 100)
x = np.linspace(0.01,.5, 10000)
xf = np.linspace(0.01, .5, 10000)
#
#F_nom = MSM_freq(H0, 5, 3, 1)
#plt.plot(H0/(4*np.pi*M), F_nom[:,0], H0/(4*np.pi*M), F_nom[:,1])
#plt.grid(True)
#plt.show()       
    

y_val = np.empty(len(x))
y_val_f = np.empty(len(xf))

#for i, j in enumerate(x):
#   y_val[i] = y(j, n, m, 610, sign)
#   
#for i, j in enumerate(xf):
#   y_val_f[i] = y_t(j, n, m, 700, sign)

#plt.figure()
#plt.plot(x, y_val)
##plt.ylim([-10,10])
#plt.grid(True)
#
#plt.figure()
#plt.plot(xf, y_val_f)
#plt.grid(True)
#plt.show()

def getF(n,s,H_0, sign):
    x = np.linspace(0.01,.5, 10000)
    y_val = np.empty(len(x))
    for i, j in enumerate(x):
        y_val[i] = y(j, n, n-s, H_0, sign)
    zeros = zero_finder(x, y_val)
    return (zeros)

F = getF(3,2, 5000, -1)
    
Happ = np.linspace(1200, 12000, 100)
Homega = Happ/(4*np.pi*M) - 1/3
res = []
for happ in Happ:
    for i, j in enumerate(x):
        temp = []
        if i == 0:
            print(happ)
        y_val_f[i] = y(j, n, n-s, happ, sign)
    zeros = zero_finder(x, y_val_f)
    res.append(zeros)
    Fres = np.array(res)
    
Freq_res = (Fres + Homega[:,np.newaxis])*4*np.pi*M*gamma    

plt.figure()   
plt.plot(Happ,Freq_res)
plt.grid(True)

plt.figure()
plt.plot(Homega, Fres)
plt.grid(True)
plt.show()


        