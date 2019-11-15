import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib

start=54562
stop=54972

def correct(data):
    import scipy
    import scipy.special
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
    

raw1=np.loadtxt('mrk421-longterm.txt',usecols=(3,1,4,5))
raw1[:,2:4]=raw1[:,2:4]*8.84
raw1=raw1[np.where(raw1[:,0]>start)]
data1=correct(raw1)
data1[:,0]=data1[:,0]-start
itp1=np.loadtxt('itplongterm.txt')
itp1[:,1]=itp1[:,1]*8.84


raw2=np.loadtxt('data/ic40-mrk501',usecols=(0,1,3,4))
raw2[:,2:4]=raw2[:,2:4]*10**7
data2=correct(raw2)
data2[:,0]=data2[:,0]-start
itp2=np.loadtxt('itms/itpmrk501.txt')

raw3=np.loadtxt('dataic40-1ES1218+304',usecols=(0,1,3,4))
raw3[:,2:4]=raw3[:,2:4]*10**7
data3=correct(raw3)
data3[:,0]=data3[:,0]-start
itp3=np.loadtxt('itms/itp1es1218+304.txt')

raw4=np.loadtxt('data/ic40-3C66A',usecols=(0,1,3,4))
raw4[:,2:4]=raw4[:,2:4]*10**7
data4=correct(raw4)
data4[:,0]=data4[:,0]-start
itp4=np.loadtxt('itms/itp3c66a.txt')

raw5=np.loadtxt('data/ic40-WComae',usecols=(0,1,3,4))
raw5[:,2:4]=raw5[:,2:4]*10**7
data5=correct(raw5)
data5[:,0]=data5[:,0]-start
itp5=np.loadtxt('itms/itpwcomae.txt')

raw6=np.loadtxt('data/ic40-1ES0806+524',usecols=(0,1,3,4))
raw6[:,2:4]=raw6[:,2:4]*10**7
data6=correct(raw6)
data6[:,0]=data6[:,0]-start
itp6=np.loadtxt('itms/itp1es0806+524.txt')
  
xform1 = matplotlib.ticker.ScalarFormatter(useOffset=False)  
xform2 = matplotlib.ticker.MultipleLocator(base=10)
cutoff=1.86


matplotlib.rc('font',family='Helvitica')
matplotlib.rcParams.update({'font.size': 16}) 
matplotlib.rcParams.update({'mathtext.default':'regular'})



baa=plt.figure(1)

baa.text(0.5, 0.01, r'MJD$-$54562', ha='center',fontweight='bold')
baa.text(0.01, 0.5, r'Flux $(10^{-11}\/cm^{-2} \/ s^{-1})}$', va='center', rotation='vertical',fontweight='bold')

aa=plt.subplot(3,1,1)
plt.errorbar(data1[:,0],data1[:,2],yerr=data1[:,3],color='blue',linestyle='none',marker='o',markersize=5)
plt.plot(itp1[:,0],itp1[:,1],color='green',linewidth=2)
plt.text(60,40,'Mrk 421',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp1[:,0],itp1[:,1],where=itp1[:,1]>6.2,alpha=0.3)
plt.ylim(ymin=0,ymax=70)
plt.xlim(xmin=-5,xmax=70)

ab=plt.subplot(3,1,2)
plt.errorbar(data1[:,0],data1[:,2],yerr=data1[:,3],color='blue',linestyle='none',marker='o',markersize=5)
plt.plot(itp1[:,0],itp1[:,1],color='green',linewidth=2)

plt.fill_between(itp1[:,0],itp1[:,1],where=itp1[:,1]>6.2,alpha=0.3)
plt.ylim(ymin=0,ymax=15)
plt.xlim(xmin=280,xmax=355)

ac=plt.subplot(3,1,3)
plt.errorbar(data1[:,0],data1[:,2],yerr=data1[:,3],color='blue',linestyle='none',marker='o',markersize=5)
plt.plot(itp1[:,0],itp1[:,1],color='green',linewidth=2)

plt.fill_between(itp1[:,0],itp1[:,1],where=itp1[:,1]>6.2,alpha=0.3)
plt.ylim(ymin=0,ymax=15)
plt.xlim(xmin=355,xmax=430)

aa.spines['right'].set_visible(False)
aa.spines['top'].set_visible(False)
aa.yaxis.set_ticks_position('left')
aa.xaxis.set_ticks_position('bottom')
aa.get_xaxis().get_major_formatter().set_useOffset(False)
aa.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
aa.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(20))
ab.spines['right'].set_visible(False)
ab.spines['top'].set_visible(False)
ab.yaxis.set_ticks_position('left')
ab.xaxis.set_ticks_position('bottom')
ab.get_xaxis().get_major_formatter().set_useOffset(False)
ab.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
ab.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(5))
ac.spines['right'].set_visible(False)
ac.spines['top'].set_visible(False)
ac.yaxis.set_ticks_position('left')
ac.xaxis.set_ticks_position('bottom')
ac.get_xaxis().get_major_formatter().set_useOffset(False)
ac.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
ac.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(5))

plt.show()


saa=plt.figure(2)
saa.text(0.5, 0.01, r'MJD$-$54562', ha='center',fontweight='bold')
saa.text(0.01, 0.5, r'Flux $(10^{-11}\/cm^{-2} \/ s^{-1})}$', va='center', rotation='vertical',fontweight='bold')

a2=plt.subplot2grid((4,3),(0,0),colspan=2)
plt.errorbar(data2[:,0],data2[:,2],yerr=data2[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp2[:,0], itp2[:,1],color='green',linewidth=2)
plt.text(40,8,'a: Mrk 501',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp2[:,0],itp2[:,1],where=itp2[:,1]>cutoff,alpha=0.3)
plt.ylim(ymin=0)
plt.xlim(xmin=-8,xmax=47)

a3=plt.subplot(4,3,3)
plt.errorbar(data2[:,0],data2[:,2],yerr=data2[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp2[:,0], itp2[:,1],color='green',linewidth=2)
plt.text(364,8,'b: Mrk 501',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp2[:,0],itp2[:,1],where=itp2[:,1]>cutoff,alpha=0.3)
plt.ylim(ymin=0)
plt.xlim(xmin=343,xmax=365)
5
a4=plt.subplot(4,1,2)
plt.errorbar(data3[:,0],data3[:,2],yerr=data3[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp3[:,0], itp3[:,1],color='green',linewidth=2)
plt.fill_between(itp3[:,0],itp3[:,1],where=itp3[:,1]>cutoff,alpha=0.3)
plt.text(290,4.7,'c: 1ES 1218+304',horizontalalignment='right',fontweight='bold')
plt.ylim(ymin=0)
plt.xlim(xmin=265,xmax=350)

a5=plt.subplot(4,1,3)
plt.errorbar(data4[:,0],data4[:,2],yerr=data4[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp4[:,0], itp4[:,1],color='green',linewidth=2)
plt.text(245,4,'d: 3C66A',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp4[:,0],itp4[:,1],where=itp4[:,1]>cutoff,alpha=0.3)
plt.ylim(ymin=0,ymax=6)
plt.xlim(xmin=161,xmax=246)

a6=plt.subplot(4,3,10)
plt.errorbar(data5[:,0],data5[:,2],yerr=data5[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp5[:,0], itp5[:,1],color='green',linewidth=2)
plt.text(16,4.5,'e: W Comae',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp5[:,0],itp5[:,1],where=itp5[:,1]>cutoff,alpha=0.3)
plt.ylim(ymin=0)
plt.xlim(xmin=-8,xmax=17)

a7=plt.subplot(4,3,11)
plt.errorbar(data5[:,0],data5[:,2],yerr=data5[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp5[:,0], itp5[:,1],color='green',linewidth=2)
plt.text(81,4.5,'f: W Comae',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp5[:,0],itp5[:,1],where=itp5[:,1]>cutoff,alpha=0.3)
plt.ylim(ymin=0)
plt.xlim(xmin=57,xmax=82)

a9=plt.subplot(4,3,12)
plt.errorbar(data6[:,0],data6[:,2],yerr=data6[:,3],color='blue',linestyle='none', marker='o',markersize=5)
plt.plot(itp6[:,0], itp6[:,1],color='green',linewidth=2)
plt.text(19,4.5,'g: 1ES 0806+524',horizontalalignment='right',fontweight='bold')
plt.fill_between(itp6[:,0],itp6[:,1],where=itp6[:,1]>cutoff,alpha=0.3)
plt.ylim(ymin=0,ymax=6)
plt.xlim(xmin=-5,xmax=20)


a2.spines['right'].set_visible(False)
a2.spines['top'].set_visible(False)
a2.yaxis.set_ticks_position('left')
a2.xaxis.set_ticks_position('bottom')
a2.get_xaxis().get_major_formatter().set_useOffset(False)
a2.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(3))
a3.spines['right'].set_visible(False)
a3.spines['top'].set_visible(False)
a3.yaxis.set_ticks_position('left')
a3.xaxis.set_ticks_position('bottom')
a3.get_xaxis().get_major_formatter().set_useOffset(False)
a3.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
a3.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(3))
a4.spines['right'].set_visible(False)
a4.spines['top'].set_visible(False)
a4.yaxis.set_ticks_position('left')
a4.xaxis.set_ticks_position('bottom')
a4.xaxis.set_major_formatter(xform1)
a4.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
a4.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(2))
a5.spines['right'].set_visible(False)
a5.spines['top'].set_visible(False)
a5.yaxis.set_ticks_position('left')
a5.xaxis.set_ticks_position('bottom')
a5.get_xaxis().get_major_formatter().set_useOffset(False)
a5.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
a5.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(2))
a6.spines['right'].set_visible(False)
a6.spines['top'].set_visible(False)
a6.yaxis.set_ticks_position('left')
a6.xaxis.set_ticks_position('bottom')
a6.get_xaxis().get_major_formatter().set_useOffset(False)
a6.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
a6.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(2))
a7.spines['right'].set_visible(False)
a7.spines['top'].set_visible(False)
a7.yaxis.set_ticks_position('left')
a7.xaxis.set_ticks_position('bottom')
a7.get_xaxis().get_major_formatter().set_useOffset(False)
a7.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
a7.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(2))
a9.spines['right'].set_visible(False)
a9.spines['top'].set_visible(False)
a9.yaxis.set_ticks_position('left')
a9.xaxis.set_ticks_position('bottom')
a9.get_xaxis().get_major_formatter().set_useOffset(False)
a9.xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(10))
a9.yaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(2))

plt.show()


int1=itp1[np.where(itp1[:,1]>6.2)]
int2=itp2[np.where(itp2[:,1]>cutoff)]
int3=itp3[np.where(itp3[:,1]>cutoff)]
int4=itp4[np.where(itp4[:,1]>cutoff)]
int5=itp5[np.where(itp5[:,1]>cutoff)]
int6=itp6[np.where(itp6[:,1]>cutoff)]

mrk421int=sum(int1[:,1]*60*60*10**-7)
mrk501int=sum(int2[:,1]*60*60*10**-7)
es1218int=sum(int3[:,1]*60*60*10**-7)
c66aint=sum(int4[:,1]*60*60*10**-7)
wcomaeint=sum(int5[:,1]*60*60*10**-7)
es0806int=sum(int6[:,1]*60*60*10**-7)


print mrk421int
otot=mrk501int+es1218int+c66aint+wcomaeint+es0806int
print otot






