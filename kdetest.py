import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats
import matplotlib

start=54562
stop=54972

def correct(data):
    #this will take negative flux measurements and correct them to be positive
    import scipy.special
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


raw1=np.loadtxt('data/ic40-1ES0806+524',usecols=(0,1,3,4))
raw1[:,2:4]=raw1[:,2:4]*10**7
data1=correct(raw1)
flux1=data1[:,2]
flux1=flux1[np.where(flux1<10)]
#itp1=np.loadtxt('itp1es0806+524.txt')

raw2=np.loadtxt('data/ic40-1ES1218+304',usecols=(0,1,3,4))
raw2[:,2:4]=raw2[:,2:4]*10**7
data2=correct(raw2)
flux2=data2[:,2]
flux2=flux2[np.where(flux2<10)]
#itp2=np.loadtxt('itp1es1218+304.txt')

raw3=np.loadtxt('data/ic40-3C66A',usecols=(0,1,3,4))
raw3[:,2:4]=raw3[:,2:4]*10**7
data3=correct(raw3)
flux3=data3[:,2]
flux3=flux3[np.where(flux3<10)]
#itp3=np.loadtxt('itp3c66a.txt')


raw6=np.loadtxt('data/ic40-WComae',usecols=(0,1,3,4))
raw6[:,2:4]=raw6[:,2:4]*10**7
data6=correct(raw6)
flux6=data6[:,2]
flux6=flux6[np.where(flux6<10)]
#itp6=np.loadtxt('itpwcomae.txt')

p1=scipy.stats.gaussian_kde(flux1)
p2=scipy.stats.gaussian_kde(flux2)
p3=scipy.stats.gaussian_kde(flux3)
p6=scipy.stats.gaussian_kde(flux6)

ax=np.linspace(-5,7,10000)


test1=p1.evaluate(ax)
test2=p2.evaluate(ax)
test3=p3.evaluate(ax)
test6=p6.evaluate(ax)


plt.figure(1)
plt.subplot(411)
plt.hist(flux1,bins=5,normed=1)
plt.plot(ax,test1)
plt.subplot(412)
plt.hist(flux2,bins=10,normed=1)
plt.plot(ax,test2)
plt.subplot(413)
plt.hist(flux3,bins=8,normed=1)
plt.plot(ax,test3)
plt.subplot(414)
plt.hist(flux6,bins=8,normed=1)
plt.plot(ax,test6)

plt.show()

def pctint(func,axis):
    m=1
    while np.trapz(func[0:m],axis[0:m])<0.5:
        m+=1  
    return axis[m]
    

print ax[np.argmax(test1)]
print ax[np.argmax(test2)]
print ax[np.argmax(test3)]
print ax[np.argmax(test6)]
avg=np.average(ax,weights=test2)
std=np.sqrt(np.average((ax-avg)**2,weights=test2))
#print avg,std


rawa=np.loadtxt('data/mrk421-longterm.txt',comments='#')
old1=rawa[np.where(rawa[:,3]>52900)]


rawb=np.loadtxt('data/mrk501_combined_lc_v0.2.slf',comments='#',usecols = (0,1,2,3,4,5,6,7,8,9,10))
rawb=rawb[0:150,:]
old2=rawb[np.where(rawb[:,3]<1)]

saa1=scipy.stats.gaussian_kde(old1[:,4])
ka1=np.linspace(-5,10,10000)
test10=saa1.evaluate(ka1)

evil=ka1[np.argmax(test10)]
print evil
 
#print ka1[m]

saa2=scipy.stats.gaussian_kde(old2[:,3])
ka2=np.linspace(-3,5,50000)
test11=saa2.evaluate(ka2)

print ka2[np.argmax(test11)]    
#print ka2[n]

matplotlib.rcParams.update({'font.size': 14})
plt.figure(2)
#plt.subplot(211)
plt.hist(old1[:,4],bins=20,normed=1)
plt.plot(ka1,test10,linewidth=2,color='black')
plt.plot((evil,evil),(0,max(test10)),color='red',linewidth=3)
plt.xlim(xmin=-1)
plt.xlabel('Blazar flux (crab units)')

#plt.subplot(212)
#plt.hist(old2[:,3],bins=20,normed=1)
#plt.plot(ka2,test11)

plt.show()










