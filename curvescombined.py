import numpy as np
import math
import matplotlib.pyplot as plt

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
            sdev=math.sqrt(((pct-z)**2+sig**2)*(1-zero)/2+sig*(z-2*pct)*math.exp(-z**2/(2*sig**2))/math.sqrt(2*math.pi))
            base[n,3]=sdev
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

def timepairs(filt,ra,dec):
    filt=filt[np.argsort(filt[:,0])]
    pairs=np.zeros((2,2))
    pairs[0,0]=ra
    pairs[0,1]=dec
    times=filt[:,0]
    pairs[1,0]=times[0]
    n=0
    m=1
    while n<len(times):
        if n+1<len(times):
            if (times[n+1]-times[n])>1.5/float(24):
                pairs[m,1]=times[n]
                pairs=np.append(pairs,np.zeros((1,2)),0)
                m+=1
                pairs[m,0]=times[n+1]
        else:
            pairs[m,1]=times[n]
        n+=1
    return pairs
    

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
#print (max(itp3[:,0])-min(itp3[:,0]))/len(itp3[:,0])*24

raw4=np.loadtxt('data/mrk421-longterm.txt',usecols=(3,1,4,5))
raw4[:,0]=raw4[:,0]-start
raw4=raw4[np.where(raw4[:,0]>-100)]
data4=correct(raw4)
data4[:,1]=data4[:,1]*8.8
itp4=np.loadtxt('itms/itplongterm.txt')
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


joined=np.concatenate((itp1,itp2,itp3,itp5,itp6),0)

fluxsort=joined[np.argsort(joined[:,1])]
fluxsort=fluxsort[::-1]
fluence=fluxsort[:,1]/24  #each time point represents an hour of observation
fluencetotal=np.sum(fluence)
curve=np.ones((len(fluence),5))
n=0
temp=0
while n<len(fluence):
    curve[n,0]=(n+1)/float(24)
    temp+=fluence[n]
    curve[n,1]=temp/fluencetotal
    curve[n,2]=fluxsort[n,0]+start
    curve[n,3]=fluxsort[n,1]
    curve[n,4]=fluxsort[n,2]
    n+=1
#that loop creates a vector with the following columns
# 0: cumulative time (days),     1: fraction of fluence
#2: mjd,     3: flux,     4: flag to separate different observations
plt.figure(1)
plt.plot(curve[:,0],curve[:,1],color='green',marker='o',markersize=2,linestyle='none')
plt.fill_between(curve[0:1261,0],curve[0:1261,1],alpha=0.3)
plt.plot([0,max(curve[:,0])],[0,1],linewidth=2,color='green')

plt.show()
curve=curve[0:15000,:]


ic40area=3.6278*(180/np.pi)**2 #area above 25 degrees
rblazar=2.3  # radius of acceptance for each blazar
ablazar=math.pi*rblazar**2 # area of all blazar in our search

ic40cnt=6668 # total neutrinos detected by ic40 above 25 degrees
nutot=ablazar*ic40cnt/ic40area # number of neutrinos expected to fall within a single blazar
nurate=nutot/(stop-start) # number of nus expected within a blazar during a single day

blzdays=curve[:,0]
ffrac=curve[:,1]

nulim=np.zeros((len(blzdays),1))
nuratio=np.zeros((len(blzdays),1))

n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.885)
    nuratio[n]=ffrac[n]/nulim[n]
    n+=1

plt.figure(2)
plt.plot(blzdays,nuratio)
plt.show()

cutoff=np.argmax(nuratio)

set1=np.zeros((1,2))
set2=np.zeros((1,2))
set3=np.zeros((1,2))
set5=np.zeros((1,2))
set6=np.zeros((1,2))
n=0
while n<cutoff:
    if curve[n,4]==1:
        set1=np.append(set1,curve[n:n+1,2:4],0)
    elif curve[n,4]==2:
        set2=np.append(set2,curve[n:n+1,2:4],0)
    elif curve[n,4]==3:
        set3=np.append(set3,curve[n:n+1,2:4],0)
    elif curve[n,4]==5:
        set5=np.append(set5,curve[n:n+1,2:4],0)
    elif curve[n,4]==6:
        set6=np.append(set6,curve[n:n+1,2:4],0)
    n+=1
        
set1=set1[1:(len(set1[:,0])-1),:]
set1=set1[np.argsort(set1[:,0])]
set2=set2[1:(len(set2[:,0])-1),:]
set2=set2[np.argsort(set2[:,0])]
set3=set3[1:(len(set3[:,0])-1),:]
set3=set3[np.argsort(set3[:,0])]
set5=set5[1:(len(set5[:,0])-1),:]
set5=set5[np.argsort(set5[:,0])]
set6=set6[1:(len(set6[:,0])-1),:]
set6=set6[np.argsort(set6[:,0])]

plt.figure(3)
plt.subplot(511)
plt.plot(set1[:,0],set1[:,1],color='blue',linestyle='none',marker='o',markersize=5)
plt.xlim(xmin=start-10,xmax=stop)
plt.ylim(ymin=0)
plt.text(stop-5,4,'1ES0806+524',horizontalalignment='right')
plt.subplot(512)
plt.plot(set2[:,0],set2[:,1],color='blue',linestyle='none', marker='o',markersize=5)
plt.xlim(xmin=start-10,xmax=stop)
plt.ylim(ymin=0)
plt.text(stop-5,4,'1ES 1218+304',horizontalalignment='right')
plt.subplot(513)
plt.plot(set3[:,0],set3[:,1],color='blue',linestyle='none', marker='o',markersize=5)
plt.xlim(xmin=start-10,xmax=stop)
plt.text(stop-5,3.5,'3C66A',horizontalalignment='right')
plt.ylim(ymin=0)
plt.subplot(514)
plt.plot(set5[:,0],set5[:,1],color='blue',linestyle='none', marker='o',markersize=5)
plt.xlim(xmin=start-10,xmax=stop)
plt.ylim(ymin=0)
plt.text(stop-5,8,'Markarian 501',horizontalalignment='right')
plt.subplot(515)
plt.plot(set6[:,0],set6[:,1],color='blue',linestyle='none', marker='o',markersize=5)
plt.xlim(xmin=start-10,xmax=stop)
plt.ylim(ymin=0)
plt.text(stop-5,4.5,'W Comae',horizontalalignment='right')

plt.show()

es0806=timepairs(set1,122.45484,52.31618)
es1218=timepairs(set2,185.34134,30.1769)
c66a=timepairs(set3,35.66505,43.03552)
mrk501=timepairs(set5,253.46757,39.76017)
wcomae=timepairs(set6,185.38204,28.23292)

print curve[cutoff,3]

#np.savetxt('1es0806+524ontimes.txt',es0806)
#np.savetxt('1es1218+304ontimes.txt',es1218)
#np.savetxt('3c66aontimes.txt',c66a)
#np.savetxt('mrk501ontimes.txt',mrk501)
#np.savetxt('wcomaeontimes.txt',wcomae)

def printer(data):
    n=0
    archive=np.zeros((len(data[:,0])-1,2),dtype='object')
    while n<len(data[:,0])-1:
        archive[n,0]='              & ' + str(round(data[n+1,0]-start,3)) + '  & '
        archive[n,1]=str(round(data[n+1,1]-start,3)) + ' \\\\'
        n+=1

    return archive
    
    

