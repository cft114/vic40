import numpy as np
import pyGPs
import math
import scipy
import matplotlib.pyplot as plt
import matplotlib

start=54562
stop=54972

def correct(data):
    #this will take negative flux measurements and correct them to be positive
    base=np.zeros((len(data[:,0]),4))
    base[:,0]=data[:,0]
    base[:,1]=data[:,1]
    base[:,2]=data[:,2]
    base[:,3]=data[:,3]
    n=0
    while n<len(data[:,2]):
        if data[n,2]<=0:
            z=data[n,2]
            sig=data[n,3]
            zero=scipy.special.erf(-z/(sig*math.sqrt(2)))
            con=1/(1-zero)
            pct=sig*math.sqrt(2)*scipy.special.erfinv(0.5/con+zero)+z
            base[n,2]=pct
        n+=1
    return base
    
def gint(data,mean=1,cov1=1,cov2=1,noise=.1):
    x=data[:,0]
    y=np.log(data[:,2])
    
    obsleng=(max(x)-min(x)+20)*24
    #z=(np.arange(0,math.ceil(obsleng),1))+min(x)
    z=(np.arange(0,math.ceil(obsleng),1)/24)+min(x)-10
  
    model = pyGPs.GPR()      # specify model (GP regression)
    model.setData(x,y)
    model.setNoise( log_sigma = np.log(noise) )
    esa=pyGPs.Core.cov.RBF(log_ell=np.log(cov1), log_sigma=np.log(cov2))
    eka=pyGPs.mean.Const(np.log(mean))
    model.setPrior(mean=eka,kernel=esa)
    model.predict(z)         # predict test cases
    yout=model.fm
    result=np.zeros((len(yout),2))
    result[:,0]=z
    result[:,1]=np.exp(yout[:,0])
    print model.nlZ
    return result

raw1=np.loadtxt('data/ic40-1ES0806+524',usecols=(0,1,3,4))
raw1[:,2:4]=raw1[:,2:4]*10**7
raw1[:,0]=raw1[:,0]-start
data1=correct(raw1)
mean1=0.6
mjd1=data1[:,0]
flux1=data1[:,2]
interpolated1=gint(data1,mean1,0.14,0.81,0.1)

raw2=np.loadtxt('data/ic40-1ES1218+304',usecols=(0,1,3,4))
raw2[:,2:4]=raw2[:,2:4]*10**7
raw2[:,0]=raw2[:,0]-start
data2=correct(raw2)
mean2=1.6
mjd2=data2[:,0]
flux2=data2[:,2]
interpolated2=gint(data2,mean2,0.55,1.25,0.08)

raw3=np.loadtxt('data/ic40-3C66A',usecols=(0,1,3,4))
raw3[:,2:4]=raw3[:,2:4]*10**7
raw3[:,0]=raw3[:,0]-start
data3=correct(raw3)
mean3=1.5
mjd3=data3[:,0]
flux3=data3[:,2]
interpolated3=gint(data3,mean3,0.89,0.86,0.19)

raw4=np.loadtxt('data/mrk421-longterm.txt',usecols=(3,1,4,5))
raw4[:,0]=raw4[:,0]-start
data4=correct(raw4)
data4[:,2:4]=data4[:,2:4]*8.84
mjd4=data4[:,0]
flux4=data4[:,2]
interpolated4=gint(data4,0.5*8.84,0.61,0.99,0.07)


raw5=np.loadtxt('data/ic40-mrk501',usecols=(0,1,3,4))
raw5[:,2:4]=raw5[:,2:4]*10**7
raw5[:,0]=raw5[:,0]-start
data5=correct(raw5)
mean5=1.5
mjd5=data5[:,0]
flux5=data5[:,2]
interpolated5=gint(data5,mean5,0.52,1.6,0.1)

raw6=np.loadtxt('data/ic40-WComae',usecols=(0,1,3,4))
raw6[:,2:4]=raw6[:,2:4]*10**7
raw6[:,0]=raw6[:,0]-start
data6=correct(raw6)
mean6=0.6
mjd6=data6[:,0]
flux6=data6[:,2]
interpolated6=gint(data6,mean6,0.99,1.63,.1)


matplotlib.rc('font',family='Helvetica')
matplotlib.rcParams.update({'font.size': 16}) 

saa=plt.figure(1)
saa.text(0.5, 0.01, r'MJD$-$54562', ha='center',fontweight='bold')
saa.text(0.01, 0.5, r'Flux $(10^{-11}\/cm^{-2} \/ s^{-1})}$', va='center', rotation='vertical',fontweight='bold')

a1=plt.subplot(717)
plt.errorbar(mjd1,flux1,yerr=data1[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(interpolated1[:,0], interpolated1[:,1],color='green',linewidth=2)
plt.text(100,3,'f. 1ES 0806+524',fontweight='bold')
plt.ylim(ymin=0)
plt.xlim(xmin=-10,xmax=stop-start)
a1.spines['right'].set_visible(False)
a1.spines['top'].set_visible(False)
a1.yaxis.set_ticks_position('left')
a1.xaxis.set_ticks_position('bottom')
a1.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(2))

a2=plt.subplot(714)
plt.errorbar(mjd2,flux2,yerr=data2[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(interpolated2[:,0], interpolated2[:,1],color='green',linewidth=2)
plt.text(100,4,'c. 1ES 1218+304',fontweight='bold')
plt.ylim(ymin=0)
plt.xlim(xmin=-10,xmax=stop-start)
a2.spines['right'].set_visible(False)
a2.spines['top'].set_visible(False)
a2.yaxis.set_ticks_position('left')
a2.xaxis.set_ticks_position('bottom')
a2.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(3))

a3=plt.subplot(715)
plt.errorbar(mjd3,flux3,yerr=data3[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(interpolated3[:,0], interpolated3[:,1],color='green',linewidth=2)
plt.text(100,3.5,'d. 3C66A',fontweight='bold')
plt.ylim(ymin=0,ymax=6)
plt.xlim(xmin=-10,xmax=stop-start)
a3.spines['right'].set_visible(False)
a3.spines['top'].set_visible(False)
a3.yaxis.set_ticks_position('left')
a3.xaxis.set_ticks_position('bottom')
a3.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(3))

a4=plt.subplot2grid((7,1),(0,0),rowspan=2)
plt.errorbar(mjd4,flux4,yerr=data4[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(interpolated4[:,0], interpolated4[:,1],color='green',linewidth=2)
plt.text(100,50,'a. Mrk 421',fontweight='bold')
plt.ylim(ymin=0,ymax=70)
plt.xlim(xmin=-10,xmax=stop-start)
a4.spines['right'].set_visible(False)
a4.spines['top'].set_visible(False)
a4.yaxis.set_ticks_position('left')
a4.xaxis.set_ticks_position('bottom')
a4.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(25))

a5=plt.subplot(713)
plt.errorbar(mjd5,flux5,yerr=data5[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(interpolated5[:,0], interpolated5[:,1],color='green',linewidth=2)
plt.text(100,9,'b. Mrk 501',fontweight='bold')
plt.ylim(ymin=0)
plt.xlim(xmin=-10,xmax=stop-start)
a5.spines['right'].set_visible(False)
a5.spines['top'].set_visible(False)
a5.yaxis.set_ticks_position('left')
a5.xaxis.set_ticks_position('bottom')
a5.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(4))

a6=plt.subplot(716)
plt.errorbar(mjd6,flux6,yerr=data6[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(interpolated6[:,0], interpolated6[:,1],color='green',linewidth=2)
plt.text(100,4.5,'e. W Comae',fontweight='bold')
plt.ylim(ymin=0)
plt.xlim(xmin=-10,xmax=stop-start)
a6.spines['right'].set_visible(False)
a6.spines['top'].set_visible(False)
a6.yaxis.set_ticks_position('left')
a6.xaxis.set_ticks_position('bottom')
a6.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(3))


plt.show()

np.savetxt('itms/itp1es0806+524.txt',interpolated1)
np.savetxt('itms/itp1es1218+304.txt',interpolated2)
np.savetxt('itms/itp3c66a.txt',interpolated3)
np.savetxt('itms/itpmrk501.txt',interpolated5)
np.savetxt('itms/itpwcomae.txt',interpolated6)
np.savetxt('itms/itpmrk421.txt',interpolated4)
