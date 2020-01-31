import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lpmn

#Input constants
a = 1.0 #sphere radius, results are independent of this value; just used for plotting
M = 1780.0 #Saturation magnetization, Oe
gamma = 2.8e6 #Gyromagnetic ratio, Hz / Gauss

#Calculate frequency
def getfreq(n,m,H_0):
    fr = np.linspace(0,100,16000)#GHz
    yr, yr2, yr3 = [], [], []
    for f in fr:
        Omega = f*1e9/(4*np.pi*gamma*M)
        H_i = H_0 - (4*np.pi/3)*M
        Omega_H = H_i/(4*np.pi*M)
        nu = Omega/(Omega_H**2-Omega**2)
        kappa = Omega_H/(Omega_H**2-Omega**2)
        xi=np.sqrt(1+1/kappa)
        yr.append(xi*lpmn(m,n,xi)[1][m][n]/lpmn(m,n,xi)[0][m][n])
        yr2.append(-1*(n+1)+m*nu)
        yr3.append(-1*(n+1)-m*nu)
    
    f = -1
    i=0
    if (yr[0]>yr2[0]):
        while(yr[i]>yr2[i]):
            i+=1
            if i>=len(yr2)-1:
                break
        if (yr[i]<yr2[i]):
            f=fr[i]
    else:
        while(yr[i]<yr2[i]):
            i+=1
            if i>=len(yr2)-1:
                break
        if (yr[i]>yr2[i]):
            f=fr[i]
    i = 0
    if (yr[0]>yr3[0]):
        while(yr[i]>yr3[i]):
            i+=1
            if i>=len(yr3)-1:
                break
        if (yr[i]<yr3[i]):
            f=fr[i]
    else:
        while(yr[i]<yr3[i]):
            i+=1
            if i>=len(yr3)-1:
                break
        if (yr[i]>yr3[i]):
            f=fr[i]
    
    print(f)
    
    if (f<0):
        print("No resonance below 100GHz")
    return f

def getfreq_2(n,m,H_0):
    fr = np.linspace(0,100,16000)#GHz
    yr, yr2, yr3 = [], [], []
    for f in fr:
        Omega = f*1e9/(4*np.pi*gamma*M)
        H_i = H_0 - (4*np.pi/3)*M
        Omega_H = H_i/(4*np.pi*M)
        nu = Omega/(Omega_H**2-Omega**2)
        kappa = Omega_H/(Omega_H**2-Omega**2)
        xi=np.sqrt(1+1/kappa)
        yr.append(xi*lpmn(m,n,xi)[1][m][n]/lpmn(m,n,xi)[0][m][n])
        yr2.append(-1*(n+1)+m*nu)
        yr3.append(-1*(n+1)-m*nu)
    
    

#Calculated coefficients and expression for magnetization found in refs [1] and [2]
def Z(m,n,xi0,nu):
    return 1.0/((n+1+xi0*lpmn(m,n,xi0)[1][m][n]/lpmn(m,n,xi0)[0][m][n])**2-(m**2)*(nu**2))

def H(m,n,xi0,nu):
    return ((a**n)*Z(m,n,xi0,nu)*(2*n+1)/lpmn(m,n,xi0)[0][m][n])*(-1.0*nu*m+(n+1+xi0*lpmn(m,n,xi0)[1][m][n]/lpmn(m,n,xi0)[0][m][n]))

def G(m,n,xi0,nu):
    return ((a**n)*Z(m,n,xi0,nu)*(2*n+1)/lpmn(m,n,xi0)[0][m][n])*((n+1+xi0*lpmn(m,n,xi0)[1][m][n]/lpmn(m,n,xi0)[0][m][n])-nu*m)

def psi(m,n,x,y,z,f,H_0):
    
    #Calculate normalized frequency and other auxiliary variables used in [1] and [2]
    Omega = f*1e9/(4*np.pi*gamma*M)
    H_i = H_0 - (4*np.pi/3)*M
    Omega_H = H_i/(4*np.pi*M)
    nu = Omega/(Omega_H**2-Omega**2)
    kappa = Omega_H/(Omega_H**2-Omega**2)
    xi0 = np.sqrt(1+1/kappa)
    
    #Convert from Cartesian to Ellipsoidal Coordinates.
    #Conversion from ellipsoidal -> Cartesian coordinates are found in reference [1], beginning of section 4
    #This is just an inversion of that set of equations
    phi = np.arctan2(y,x)
    R2 = (1.0+kappa)*(z**2)+(x**2)+(y**2)+kappa*(a**2)
    xi = np.sqrt( ( R2 + np.sqrt((R2**2)-4*(a**2)*kappa*(1.0+kappa)*(z**2)) )/( 2*(a**2)*kappa ) )*np.sign(z)
    
    #ref [2] eqn 58
    ceta = (z/(a*xi))*np.sqrt((1.0+kappa)/kappa)
    
    return lpmn(m,n,xi)[0][m][n]*lpmn(m,n,ceta)[0][m][n]*(G(m,n,xi0,nu)*np.cos(m*phi)+1j*H(m,n,xi0,nu)*np.sin(m*phi))

def dpsidx(m,n,x,y,z,f,H_0):
    return (psi(m,n,x+0.001,y,z,f,H_0)-psi(m,n,x,y,z,f,H_0))/0.001

def dpsidy(m,n,x,y,z,f,H_0):
    return (psi(m,n,x,y+0.001,z,f,H_0)-psi(m,n,x,y,z,f,H_0))/0.001

def mx(m,n,x,y,z,f,H_0):
    Omega = f*1e9/(4*np.pi*gamma*M)
    H_i = H_0 - (4*np.pi/3)*M
    Omega_H = H_i/(4*np.pi*M)
    nu = Omega/(Omega_H**2-Omega**2)
    kappa = Omega_H/(Omega_H**2-Omega**2)
    return (kappa*dpsidx(m,n,x,y,z,f,H_0)-nu*dpsidy(m,n,x,y,z,f,H_0))

def my(m,n,x,y,z,f,H_0):
    Omega = f*1e9/(4*np.pi*gamma*M)
    H_i = H_0 - (4*np.pi/3)*M
    Omega_H = H_i/(4*np.pi*M)
    nu = Omega/(Omega_H**2-Omega**2)
    kappa = Omega_H/(Omega_H**2-Omega**2)
    return (1j*nu*dpsidx(m,n,x,y,z,f,H_0)+kappa*dpsidy(m,n,x,y,z,f,H_0))


#Plot the frequency of first 5 (n,n,0) modes
#Brange=np.linspace(0,8000,10)
#freqs=[]
#ns=[1,2,3,4,5]
#for n in ns:
#    freqsi=[]
#    for B in Brange:
#        freqsi.append(getfreq(n,n,B))
#    freqs.append(freqsi)
#plt.figure(1)
#make nice plot
#matplotlib.rcParams['axes.linewidth']=1.5
#matplotlib.rcParams.update({'figure.autolayout': True})
#plt.tick_params(which='both', direction='in', length=12, width=1.5, right=True, top=True, labelsize=22)
#plt.tick_params(which='minor',length=5)
#plt.tick_params(axis='y',pad=10)
#plt.tick_params(axis='x',pad=10)
#plt.minorticks_on()
#for fr in range(len(freqs)):
#    plt.plot(Brange,freqs[fr],label="(%i,%i,0)"%(ns[fr],ns[fr]),linewidth=2)
#plt.ylim([0,20])
#plt.ylabel("Mode Frequency (GHz)",fontsize=20)
#plt.xlim([0,4000])
#plt.xlabel("External Field (Gauss)",fontsize=20)
#plt.legend(fontsize=16,frameon=False,loc='upper left')
#plt.show()


#Plot the magnetization pattern of a given mode

def plot_mode(n,m,H0):
    mxs, mys = [], []
    zs = np.linspace(-0.80,0.80,12)
    zs[11] = 0
    print(zs)
    xs=[]
    ys=[]
    f=getfreq(n,m,H0)
    for z in zs:
        mxcur, mycur  = [], []
        xsc, ysc = [], []
        for r in np.linspace(0.1,np.sqrt(a**2-z**2)-0.1,5):
            for theta in np.linspace(0.08,2*np.pi+0.08,3+int((r/a)*36)):
                x = r*np.cos(theta)
                y = r*np.sin(theta)
                if (x**2+y**2+z**2<a**2):
                    xsc.append(x)
                    ysc.append(y)
                    mxcur.append(np.real(mx(m,n,x,y,z,f,H0)))
                    mycur.append(np.real(my(m,n,x,y,z,f,H0)))
        mxs.append(mxcur)
        mys.append(mycur)
        xs.append(xsc)
        ys.append(ysc)
    
    fig, axs = plt.subplots(3,4)
    fig.suptitle("Mode (n,m,r) = (%s,%s,0)"%(n,m),fontsize=24)
    xplot=0
    yplot=0
    for i in range(len(zs)):
        axs[yplot,xplot].quiver(xs[i],ys[i],mxs[i],mys[i],units="xy",facecolor='white',edgecolor='black',lw=1,pivot='middle',headwidth=10,headlength=5,minshaft=8)
        axs[yplot,xplot].set_aspect(1.0)
        axs[yplot,xplot].set_xlim([-1.0,1.0])
        axs[yplot,xplot].set_ylim([-1.0,1.0])
        axs[yplot,xplot].add_artist(plt.Circle((0, 0), np.sqrt(a**2-zs[i]**2), fill=False, color='black',linewidth=3))
        xplot+=1
        if xplot>=4:
            xplot=0
            yplot+=1
    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    plt.gcf().subplots_adjust(top=0.55)
    plt.tight_layout()
    plt.show()

plot_mode(2,2,2000)

#References:
#[1] Roschmann and Dotsch, phys. stat. sol. (b) 82, 11 (1977); https://doi.org/10.1002/pssb.2220820102
#[2] Fletcher and Bell, Journal of Applied Physics 30, 687 (1959); https://doi.org/10.1063/1.1735216