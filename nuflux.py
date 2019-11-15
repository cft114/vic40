import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

    
def getc(phi,gamma,th,nugam):
    c=phi*(1-gamma)/(30**(1-gamma)-th**(1-gamma))
    spiral=c/(4*(0.5**(2-nugam)))
    return spiral
    

#threshold energies
t1=0.4
t2=0.3
t3=0.2
t4=0.2
t5=0.2
t6=0.3  
      
#nu spectral indices
g1=2.0
g2=2.5
g3=2.2
g4=2.0
g5=3.0
g6=2.7
          
c1=getc(4.17,2.2,t1,g1) # mrk 421

c2=getc(0.602,2.72,t2,g2) # mrk 501

c3=getc(0.308,3.1,t3,g3) # 1es 1218+304


c4=getc(0.25,4.1,t4,g4) # 3c 66a

c5=getc(0.109,3.8,t5,g5) #wcomae

c6=getc(0.00209,3.6,t6,g6) #1es0806+524



npoints=1000000
cap=1000
z=np.arange(0,npoints+1,1)*cap/float(npoints)
dx=z[1]-z[0]

aea=np.loadtxt('ic40area.txt')
en=aea[:,0]
ar=aea[:,1]

itp=interp1d(en, ar, kind='cubic')
dat=itp(z[np.where(z>0.2)])
plt.loglog(z[np.where(z>0.2)],dat)
plt.loglog(en,ar,linestyle='none',marker='x')
plt.show()


es0806cnu=9.5*(10**-4)


x1=z[np.where(z>t1)]
mrk421nu=np.trapz(c1*itp(x1)*x1**(-g1),dx=dx)

x2=z[np.where(z>t2)]
mrk501nu=np.trapz(c2*itp(x2)*x2**(-g2),dx=dx)

x3=z[np.where(z>t3)]
es1218nu=np.trapz(c3*itp(x3)*x3**(-g3),dx=dx)

x4=z[np.where(z>t4)]
c66anu=np.trapz(c4*itp(x4)*x4**(-g4),dx=dx)

x5=z[np.where(z>t5)]
wcomaenu=np.trapz(c5*itp(x5)*x5**(-g5),dx=dx)

x6=z[np.where(z>t6)]
es0806nu=np.trapz(c6*itp(x6)*x6**(-g6),dx=dx)



nusum=es0806nu+es1218nu+c66anu+mrk501nu+wcomaenu





