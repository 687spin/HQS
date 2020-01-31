# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:09:04 2016

@author: iwh
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
import scipy.constants as sc
import math


# Define Constants for the simulation
w = 16e-6   # CPW Center Conductor Width 
h = 7e-7    # CPW Height
gap = 9e-6  # CPW Gap Distance
gnd = 25e-6 # CPW Ground Plane Width
theta = 15/180 * np.pi  # Read Head Angle of Attack
current = .05           # Current in CPW center conductor
int_resolution = 100    # Resolution for Numerical integration if using np.arragne. If Using real data, resolution set by import array size
p_depth = 40e-9         # London Penetration Depth
fh = 1.0e-6
x_scale = 1


#j = current / (w * h)
beta = sc.mu_0/(2*sc.pi)
#gamma = beta * j(a,superconducting)
superconducting = True


def j(x, superconducting):
    if superconducting:
       return (current / (w * h))*(1 - ((2 * x)/w)**2)**(-0.5)
    else:    
       return  current / (w * h)
        
def gamma(x, superconducting):
    return beta * j(x, superconducting)


readpath = r'X:\Haygood\HD_RH\16um Nb 4mm CPW\20160629_AC_47mA_10KHz_10K_Step1_flatter.csv'
data = np.genfromtxt(readpath, delimiter = ',', skip_header = 1)

def ctr_cnd(x, w, h):
    res = np.zeros_like(x)
    for i, val in enumerate(x):
        if abs(val) <= w/2:
            res[i] = h*100
        else:
            res[i] = 0
    return res

def cpw_str(x, w, h, gap, gnd):
    res = np.zeros_like(x)
    for i, val in enumerate(x):
        if abs(val) > w/2 + gap + gnd:
            res[i] = 0
        elif abs(val) <=w/2 + gap + gnd and abs(val) >= w/2 + gap:
            res[i] = h * 50
        elif abs(val) >= w/2 and abs(val) < w/2 + gap:
            res[i] = 0
        else:
            res[i] = h * 50
    return res
    


def dhx(z, x, xo, zo, superconducting):
    return gamma(x, superconducting) * ((zo + z)/ ((xo - x)**2 + (zo + z )**2))
    
def dhz(z, x, xo, zo, superconducting):
    return gamma(x, superconducting) * ((x - xo) / ((xo - x)**2 + (zo + z )**2))
    
def H(xo, zo):
    return integrate.quad(dhx)


    
    
#ans, err = integrate.quad(simple1, 0, sc.pi)
#print (ans, err)

a = np.arange(-4*w, 4*w, w/int_resolution)

'''
def int_simple1(a):
    res = np.zeros_like(a)
    for i, val in enumerate(a):
        y, err = integrate.quad(simple1, 0, val, )
        res[i] = y
    return res
'''

    
def int_simple1(a):
    res = np.zeros_like(a)
    for i, val in enumerate(a):
        y, err = integrate.quad(simple1, 0, val, args = (val,) )
        res[i] = y
    return res


def hx(a, w, h, fh, superconducting):
    res = np.zeros_like(a)
    for i, val in enumerate(a):
        y, err = integrate.dblquad(dhx, -w/2, w/2, lambda x: 0, lambda x: h, args = (val, fh, superconducting) )
        res[i] = y 
    return res
    
def hz(a, w, h, fh, superconducting):
    res = np.zeros_like(a)
    for i, val in enumerate(a):
        y, err = integrate.dblquad(dhz, -w/2, w/2, lambda x: 0, lambda x: h, args = (val, fh, superconducting) )
        res[i] = y 
    return res
    

def sim_sig(x, w, h, superconducting, fh, theta, x_scale):
    x = x * x_scale
    return hx(x, w, h, fh, superconducting) * np.sin(theta) + hz(x, w, h, fh, superconducting) * np.cos(theta)
 
def residuals(p, w, h, superconducting, x, y):
    fh, theta, x_scale = p
    err = y - sim_sig(x, w, h, superconducting, fh, theta, x_scale)
    return err        
    
sim_sig_N = sim_sig(a, w, h, False, fh, theta, x_scale)
sim_sig_SC = sim_sig(a, w, h, True, fh, theta, x_scale)     
#sim_sig_SC = hx(a, w, h, fh, superconducting) * np.sin(theta) + hz(a, w, h, fh, superconducting) * np.cos(theta)
#sim_sig_N = hx(a, w, h, fh, False) * np.sin(theta) + hz(a, w, h, fh, False) * np.cos(theta)



datay = data[:,1]/max(data[:,1])*max(sim_sig_N)
deltax_sim = -(np.argmax(sim_sig_N) - np.argmin(sim_sig_N))
deltax_raw = -(np.argmax(datay) - np.argmin(datay))
datax = (data[:,0] + len(data[:,0])/2) *(w / deltax_raw) 
deltax_pos = -(datax[np.argmax(datay)] - np.argmax(sim_sig_N) * w/(2 * int_resolution))/2
datax = ((data[:,0] + len(data[:,0])/2) *(w / deltax_raw)) + deltax_pos
sim_sig_N = sim_sig(datax, w, h, False, fh, theta, x_scale)
sim_sig_SC = sim_sig(datax, w, h, True, fh, theta, x_scale) 
#plt.plot(a, hx(a, w, h, superconducting), a, hz(a, w, h, superconducting), a, ctr_cnd(a, w, h))

p0 = [1e-6, 14*np.pi/180, 1.0]
np.array(p0)
plsq = optimize.leastsq(residuals, p0, args = (w, h, False, datax, datay))
plsq_fit = plsq[0].tolist()

sim_sig_ft = sim_sig(datax, w, h, False, plsq_fit[0], plsq_fit[1], plsq_fit[2])



plt.plot(datax, sim_sig_ft, color = 'r', linewidth = 2)
plt.plot(datax, sim_sig_SC, color = 'b', linewidth = 2)
plt.plot(datax, cpw_str(datax, w, h, gap, gnd), color = '0.75')
plt.plot(datax, datay, color = 'k', linestyle = '--')
plt.axvline(x = -w/2, color = 'k')
plt.axvline(x = w/2, color = 'k')
plt.axvline(x = 0, color = 'k')
plt.axhline(y = 0, color = 'k')
plt.xlim(-4*w, 4*w)
#plt.axvline(x = -w/2, color = 'k')
#plt.axvline(x = w/2, color = 'k')
print(deltax_pos)
