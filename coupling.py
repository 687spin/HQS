import numpy as np
from spheremodes_backup import mx, my, getfreq
import matplotlib.pyplot as plt

cellsize=1.0
a=1.0
scale = 1

A = np.genfromtxt("./results/CPWplot_A.txt")

#sphere center coordinates and radius in units used in vector potential simulation
xcenter = 260
ycenter = 203
radius = 100

plt.figure(1)
plt.title("Vector Potential, $A_{z}$")
plt.contourf(A.T,levels=35)
plt.axes().add_artist(plt.Circle((xcenter,ycenter), radius, fill=False, color='black'))
plt.xlabel("x ($\mu m$)")
plt.ylabel("y ($\mu m$)")
plt.show()

def eta(n,m,L):
    H=6000
    freq=getfreq(n,m,H)
    ms=[]#magnetization (x- and y- comp only, no z-comp in RF magnetic field so unnecessary)
    hs=[]#external magnetic field
    eta=0
    for r in np.linspace(0.01,a,20):
        for theta in np.linspace(0,np.pi,20):
            for phi in np.linspace(0,2*np.pi,20):
                #sphere coordinates
                x=r*np.sin(theta)*np.cos(phi)
                y=r*np.sin(theta)*np.sin(phi)
                z=r*np.cos(theta)
                dr=(a-0.01)/10.0
                dtheta=np.pi/10.0
                dphi=2*np.pi/20.0
                
                mag=[mx(m,n,x,y,z,freq,H)*np.cos(np.pi*z/L),my(m,n,x,y,z,freq,H)*np.cos(np.pi*z/L)]
                ms.append(mag[0]**2+mag[1]**2)
                
                #get field coordinates corresponding to sphere coordinates
                xf=xcenter+int(x*(radius/a))
                yf=ycenter+int(y*(radius/a))
                h=[(A[xf,yf+1]-A[xf,yf]),(A[xf,yf]-A[xf+1,yf])]
                hs.append(h[0]**2+h[1]**2)
                eta+=dr*dtheta*dphi*(h[0]*mag[0]+h[1]*mag[1])*r*np.sin(theta)
    eta = np.abs(3*eta/(np.sqrt(freq)*np.max(hs)*np.max(ms)*4*np.pi*(a**3)))
    return eta

etas=[eta(n,n,100) for n in range(1,5)]
plt.figure(2)
plt.scatter(range(1,5),[eta/etas[0] for eta in etas])
plt.show()
#print(etas[0][0])
#
#plt.figure(2)
#plt.title("$\\eta$")
#plt.xlabel("n")
#plt.ylabel("$\\eta$ (A.U.)")
#for i in range(len(etas)):
#    etas[i] = [etas[i][j]/scale for j in range(len(etas[i]))]
#    plt.plot([j+i for j in range(len(etas[i]))],etas[i],color="C%i"%i)
#    plt.scatter([j+i for j in range(len(etas[i]))],etas[i],color="C%i"%i,label="(n,n,0)"%i)
#plt.legend()
#plt.show()