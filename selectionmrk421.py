import numpy as np
import math
import matplotlib.pyplot as plt
import scipy

start=54562
stop=54972

def correct(data):
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

def poiprob(clim,lam):
    n=0
    ptot=0
    while n<=clim:
        ptot+=(lam**n)*math.exp(-lam)/math.factorial(n)
        n+=1
    return ptot
        

def poilim(lam,pval):
    if 0<pval<0.5:
        ctest=1+math.floor(lam)
        ptest=poiprob(ctest,lam)
        while ptest>pval:
            ctest-=1
            ptest=poiprob(ctest,lam)
        clim=ctest
    elif 0.5<=pval<1:
        ctest=math.floor(lam)
        ptest=poiprob(ctest,lam)
        while ptest<pval:
            ctest+=1
            ptest=poiprob(ctest,lam)
        clim=ctest+1
    return clim

raw1=np.loadtxt('data/ic40-1ES0806+524',usecols=(0,1,3,4))
raw1[:,2:4]=raw1[:,2:4]*10**7
data1=correct(raw1)
itp1=np.loadtxt('itms/itp1es0806+524.txt')
itp1=np.concatenate((itp1,np.ones((len(itp1[:,0]),1))),1)
#print (max(itp1[:,0])-min(itp1[:,0]))/len(itp1[:,0])*24

raw2=np.loadtxt('data/ic40-1ES1218+304',usecols=(0,1,3,4))
raw2[:,2:4]=raw2[:,2:4]*10**7
data2=correct(raw2)
itp2=np.loadtxt('itms/itp1es1218+304.txt')
itp2=np.concatenate((itp2,2*np.ones((len(itp2[:,0]),1))),1)
#print (max(itp2[:,0])-min(itp2[:,0]))/len(itp2[:,0])*24

raw3=np.loadtxt('data/ic40-3C66A',usecols=(0,1,3,4))
raw3[:,2:4]=raw3[:,2:4]*10**7
data3=correct(raw3)
itp3=np.loadtxt('itms/itp3c66a.txt')
itp3=np.concatenate((itp3,3*np.ones((len(itp3[:,0]),1))),1)
#print (max(itp3[:,0])-min(itp3[:,0]))/len(itp3[

raw4=np.loadtxt('data/mrk421-longterm.txt',usecols=(3,1,4,5))
raw4[:,0]=raw4[:,0]-start
raw4=raw4[np.where(raw4[:,0]>-100)]
data4=correct(raw4)
data4[:,1]=data4[:,1]*8.8
itp4=np.loadtxt('itms/itplongterm.txt')
itp4[:,1]=itp4[:,1]*8.8
itp4=np.concatenate((itp4,4*np.ones((len(itp4[:,0]),1))),1)

raw5=np.loadtxt('data/ic40-mrk501',usecols=(0,1,3,4))
raw5[:,2:4]=raw5[:,2:4]*10**7
data5=correct(raw5)
itp5=np.loadtxt('itms/itpmrk501.txt')
itp5=np.concatenate((itp5,5*np.ones((len(itp5[:,0]),1))),1)
#print (max(itp5[:,0])-min(itp5[:,0]))/len(itp5[:,0])*24

raw6=np.loadtxt('data/ic40-WComae',usecols=(0,1,3,4))
raw6[:,2:4]=raw6[:,2:4]*10**7
data6=correct(raw6)
itp6=np.loadtxt('itms/itpwcomae.txt')
itp6=np.concatenate((itp6,6*np.ones((len(itp6[:,0]),1))),1)
#print (max(itp6[:,0])-min(itp6[:,0]))/len(itp6[:,0])*24

joined=np.concatenate((itp1,itp2,itp3,itp4,itp5,itp6),0)


fluxsort=joined[np.argsort(joined[:,1])]
fluxsort=fluxsort[::-1]
fluence=fluxsort[:,1]/24  #each time point represents an hour of observation, to better than one part in a thousand
fluencetotal=np.sum(fluence)
curve=np.ones((len(fluence),5))
n=0
temp=0
while n<len(fluence):
    curve[n,0]=(n+1)/float(24)  #running ttime count
    temp+=fluence[n]
    curve[n,1]=temp/fluencetotal   #fraction of total fluence
    curve[n,2]=fluxsort[n,0]+start  # mjd of measurement
    curve[n,3]=fluxsort[n,1]  # flux of measurement
    curve[n,4]=fluxsort[n,2]  # sorting flag
    n+=1

plt.figure(1)
plt.plot(curve[:,0],curve[:,1],color='green',marker='o',markersize=2,linestyle='none')
plt.plot([0,1526],[0,1],color='green',linewidth=2)
#plt.fill_between(curve[0:636,0],curve[0:636,1],alpha=1)
#plt.fill_between(curve[0:6317,0],curve[0:6317,1],alpha=.3,color='red')
plt.xlim(xmin=0,xmax=1526)
plt.xlabel('time (days)')
plt.ylabel('fraction of fluence')

plt.show()
curve=curve[0:10000,:]

ic40area=3.6278*(180/np.pi)**2 #area above 25 degrees
rblazar=2.3  # radius of acceptance for each blazar
ablazar=np.pi*rblazar**2 # area of all blazar in our search

ic40cnt=6668 # total neutrinos detected by ic40
nutot=ablazar*ic40cnt/ic40area # number of neutrinos expected to fall within a single blazar
nurate=nutot/(stop-start) # number of nus expected within a blazar during a single day

blzdays=curve[:,0:1]
ffrac=curve[:,1]

nulim=np.zeros((len(blzdays),1))
nuratio=np.zeros((len(blzdays),1))

n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.995)
    nuratio[n]=ffrac[n]/(nulim[n])
    n+=1

plt.figure(2)
plt.plot(blzdays,nuratio)
plt.show()

cutoff=np.argmax(nuratio)

mrk421=np.zeros((1,2))

n=0
while n<cutoff:
    mrk421=np.append(mrk421,curve[n:n+1,2:4],0)
    n+=1
        
mrk421=mrk421[1:(len(mrk421[:,0])-1),:]
mrk421=mrk421[np.argsort(mrk421[:,0])]

plt.figure(3)
plt.plot(mrk421[:,0],mrk421[:,1],color='blue',linestyle='none', marker='o',markersize=5)
plt.xlim(xmin=start,xmax=stop)
plt.ylim(ymin=0)
plt.text(stop-5,80,'Markarian 421',horizontalalignment='right')
plt.show()

ra=166.11381
dec=38.20883
pairs=np.zeros((2,2))
pairs[0,0]=ra
pairs[0,1]=dec
times=mrk421[:,0]
pairs[1,0]=times[0]#+start
n=0
m=1
while n<len(times):
    if n+1<len(times):
        if (times[n+1]-times[n])>1.5/float(24):
            pairs[m,1]=times[n]#+start
            pairs=np.append(pairs,np.zeros((1,2)),0)
            m+=1
            pairs[m,0]=times[n+1]#+start
    else:
        pairs[m,1]=times[n]
    n+=1



