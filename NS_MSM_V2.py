# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 12:44:45 2019

@author: iwh
"""

import numpy as np
import matplotlib.pyplot as plt
import MSM_Zeros_SI

a = 1
nx = 50
ny = 50
sign = 1

nr = 5
nphi = 24
nz = 7

n =2
s = 0
r = 0
if sign == -1:
    r = r-1
M = 141 #Saturation magnetization, Oe
gamma = 2.8e6 #Gyromagnetic ratio, Hz / Gauss
H_0 = 5000

def psi(n,s,r, x, y, z, f, H_0, sign):
    #F = MSM_Zeros_SI.getF(n, s, H_0)[r]
    freq = (F + (H_0/(4*np.pi*M) - 1/3))*4*np.pi*M*gamma
    OmH = H_0/(4*np.pi*M) - (1/3)
    Omega = F + OmH
    kappa = -OmH/(F*(F+2*OmH))
    #print('psi_kappa: {}'.format(kappa))
    nu = -(F + OmH)/(F*(F+2*OmH))
    #print('psi_nu: {}'.format(nu))
    xi = np.sqrt(1 + 1/kappa)
    m = (n - s)*sign
    print(m)
    if s == 0:
        p_func = np.power((x + 1j*y),m)
        return p_func
    elif s == 1:
        p_func = np.power((x + 1j*y),m)*z
        return p_func
    elif s ==2:
        p_func = np.power((x + 1j*y),m)*((2*((m+1)*(2*m+3))*np.power(z,2)*xi**2) - ((2*m + 3)*(x**2 + y**2)/kappa) - (6*(m + 1)))
        return p_func
    elif s ==3:
        p_func = np.power((x + 1j*y),m)*z*((2*((m+1)*(2*m+5))*np.power(z,2)*xi**2) - 3*((2*m + 5)*(x**2 + y**2)/kappa) - (6*(m + 1)))
        return p_func
      
        
def dpsidx(n,s,r,x,y,z,f,H_0, sign):
    return (psi(n,n-s,r,x-0.0005,y,z,f,H_0, sign)-psi(n,n-s,r,x+0.0005,y,z,f,H_0, sign))/0.001

def dpsidy(n,s,r,x,y,z,f,H_0, sign):
    return (psi(n,n-s,r,x,y-0.0005,z,f,H_0, sign)-psi(n,n-s,r,x,y+0.0005,z,f,H_0, sign))/0.001

def mx(n,s,r,x,y,z,F,H_0, sign):
    OmH = H_0/(4*np.pi*M) - (1/3)
    Omega = F + OmH
    kappa = -OmH/(F*(F+2*OmH))
    nu = -(F + OmH)/(F*(F+2*OmH))
    xi = np.sqrt(1 + 1/kappa)
    return (kappa*dpsidx(n,n-s,r,x,y,z,F,H_0, sign)-1j*nu*dpsidy(n,n-s,r,x,y,z,F,H_0, sign))

def my(n,s,r,x,y,z,F,H_0, sign):
    OmH = H_0/(4*np.pi*M) - (1/3)
    Omega = F + OmH
    kappa = -OmH/(F*(F+2*OmH))
    nu = -(F + OmH)/(F*(F+2*OmH))
    xi = np.sqrt(1 + 1/kappa)
    return (1j*nu*dpsidx(n,n-s,r,x,y,z,F,H_0, sign)+kappa*dpsidy(n,n-s,r,x,y,z,F,H_0, sign))

#sphere_points_xy = np.empty((nr, ntheta, nphi))

t_step = 180/(nz+1)
theta = np.arange(t_step, t_step*(nz+1), t_step)
zp_vals = zo = a*np.cos(theta*np.pi/180)

x_vals = np.zeros((nr*nphi, nz))
y_vals = np.zeros((nr*nphi, nz))
z_vals = np.zeros((nr*nphi, nz))
mx_vals = np.zeros((nr*nphi, nz))
my_vals = np.zeros((nr*nphi, nz))
phi_vals = np.arange(0,360, 360/nphi)
plot_r_vals = a*np.sin(theta*np.pi/180)
print (MSM_Zeros_SI.getF(n, s, H_0, sign))
F = MSM_Zeros_SI.getF(n, s, H_0, sign)[r]
for i in range(len(theta)):
    ro = a*np.sin(theta[i]*np.pi/180)
    zo = a*np.cos(theta[i]*np.pi/180)
    r_vals = np.linspace(ro*0.25, ro, nr)
    for ii in range(len(r_vals)):
       for iii in range(len(phi_vals)):
           #print (ii)
           #print(iii)
           x = r_vals[ii]*np.cos(phi_vals[iii]*np.pi/180)
           y = r_vals[ii]*np.sin(phi_vals[iii]*np.pi/180)
           if zo <= 1e-6 and zo >= -1e-6:
               z = 0
           else:
               z = zo    
           x_vals[ii*nphi+iii,i] = x
           y_vals[ii*nphi+iii,i] = y
           z_vals[ii*nphi+iii,i] = z
           mx_vals[ii*nphi+iii,i] = np.real(mx(n,s,r,x,y,z,F, H_0, sign))
           my_vals[ii*nphi+iii,i] = np.real(my(n,s,r,x,y,z,F, H_0, sign))
           

#%%
           
           
           
phi_val = 90
ra = a*np.sin(phi_val*np.pi/180)
ntheta = 24
theta_vals = np.linspace(0,180, 13)
lz_vals = ra*np.cos(theta_vals*np.pi/180)
lz_vals[6] = 0 
nx_vals = []
nz_vals = []
nmx_vals = []
nmz_vals = []
phi_val = 112.5
for i, j in enumerate(lz_vals):
    y = 0
    if i == 0 or i == 12:
        nx_vals.append(0)
        x = 0
        nz_vals.append(j)
        z = j
        v_val = np.real(mx(n,s,r,x,y,z,F, H_0, sign))
        nmx_vals.append(v_val)
        nmz_vals.append(0)
    elif i == 1 or i == 11:
        sx = ra*np.sin(theta_vals[i]*np.pi/180)
        sx_vals = np.linspace(-sx, sx, 2)
        for ii in sx_vals:
            nx_vals.append(ii)
            x = ii
            nz_vals.append(j)
            z = j
            v_val = np.real(mx(n,s,r,x,y,z,F, H_0, sign))
            nmx_vals.append(v_val)
            nmz_vals.append(0)
    elif i == 2 or i == 10:
        sx = ra*np.sin(theta_vals[i]*np.pi/180)
        sx_vals = np.linspace(-sx, sx, 5)
        for ii in sx_vals:
            nx_vals.append(ii)
            x = ii
            nz_vals.append(j)
            z = j
            v_val = np.real(mx(n,s,r,x,y,z,F, H_0, sign))
            nmx_vals.append(v_val)
            nmz_vals.append(0)
    else:
        sx = ra*np.sin(theta_vals[i]*np.pi/180)
        sx_vals = np.linspace(-sx, sx, 7)
        for ii in sx_vals:
            nx_vals.append(ii)
            x = ii
            nz_vals.append(j)
            z = j
            v_val = np.real(mx(n,s,r,x,y,z,F, H_0, sign))
            nmx_vals.append(v_val)
            nmz_vals.append(0)
        
fig = plt.figure()
ax = fig.add_subplot(111)

ax.quiver(nx_vals, nz_vals, nmx_vals, nmz_vals, pivot = 'middle')
ax.set_aspect(1.0)
#ax.set_xlim([-1.2, 1.2])
#ax.set_ylim([-1.2, 1.2])
ax.set_yticklabels([])
ax.set_xticklabels([])
#ax.add_artist(plt.Circle((0, 0),2*plot_r_vals[3], fill=False, color='black',linewidth=1))
ax.add_artist(plt.Circle((0, 0),ra, fill=False, color='black',linewidth=1))
plt.text(0,0.125, 'xz', fontsize = 36, ha = 'center', va = 'center', bbox = dict(facecolor = 'none', edgecolor = 'none'))
plt.title('({},{},{})'.format(n, n-s, r), fontsize = 36)
plt.xticks([])
plt.yticks([])   
plt.show()     
       

           


#%%
 
zpl =2
           
fig = plt.figure()
ax = fig.add_subplot(111)

ax.quiver(x_vals[:,zpl], y_vals[:,zpl], mx_vals[:,zpl], my_vals[:,zpl], pivot = 'middle')
ax.set_aspect(1.0)
ax.set_xlim([-1.2, 1.2])
ax.set_ylim([-1.2, 1.2])
ax.set_yticklabels([])
ax.set_xticklabels([])
#ax.add_artist(plt.Circle((0, 0),2*plot_r_vals[3], fill=False, color='black',linewidth=1))
ax.add_artist(plt.Circle((0, 0),plot_r_vals[3], fill=False, color='black',linewidth=1))
plt.text(0, 0, 'z={}'.format(round(zp_vals[zpl],1)), fontsize = 24, ha = 'center', va = 'center')
plt.text(-1.1,1, 'xy', fontsize = 36)
plt.title('({},{},{})'.format(n, n-s, r), fontsize = 36)
plt.xticks([])
plt.yticks([])
#plt.axis('off')
   


def mxy(n,m,r, x, y):
    if n == m:
        mag_x = m*np.power(complex(x,y), m-1)
        mag_y = m*complex(0,np.power(complex(x,y),m-1))
        return mag_x, mag_y
    
    

    
