import numpy as np
import matplotlib.pyplot as plt

start=54562
stop=54972

def poiprob(clim,lam):
    n=0
    ptot=0
    while n<=clim:
        ptot+=(lam**n)*np.exp(-lam)/np.math.factorial(n)
        n+=1
    return ptot
        

def poilim(lam,pval):
    if 0<pval<0.5:
        ctest=1+np.floor(lam)
        ptest=poiprob(ctest,lam)
        while ptest>pval:
            ctest-=1
            ptest=poiprob(ctest,lam)
        clim=ctest
    elif 0.5<=pval<1:
        ctest=np.floor(lam)
        ptest=poiprob(ctest,lam)
        while ptest<pval:
            ctest+=1
            ptest=poiprob(ctest,lam)
        clim=ctest+1
    return clim
    
#import light curves, add marker to each blazar
itp1=np.loadtxt('itp1es0806+524.txt')
itp1=np.concatenate((itp1,np.ones((len(itp1[:,0]),1))),1)

itp2=np.loadtxt('itp1es1218+304.txt')
itp2=np.concatenate((itp2,2*np.ones((len(itp2[:,0]),1))),1)

itp3=np.loadtxt('itp3c66a.txt')
itp3=np.concatenate((itp3,3*np.ones((len(itp3[:,0]),1))),1)


itp4=np.loadtxt('itplongterm.txt')
itp4[:,1]=itp4[:,1]*8.8
itp4=np.concatenate((itp4,4*np.ones((len(itp4[:,0]),1))),1)

itp5=np.loadtxt('itpmrk501.txt')
itp5=np.concatenate((itp5,5*np.ones((len(itp5[:,0]),1))),1)

itp6=np.loadtxt('itpwcomae.txt')
itp6=np.concatenate((itp6,6*np.ones((len(itp6[:,0]),1))),1)

joined=np.concatenate((itp1,itp2,itp3,itp4,itp5,itp6),0)


fluxsort=joined[np.argsort(joined[:,1])] #sort by fluence
fluxsort=fluxsort[::-1] #brightets first
fluence=fluxsort[:,1]/24  #each time point represents an hour of observation, to better than one part in a thousand
fluencetotal=np.sum(fluence) #total fluence
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
plt.plot([0,1785],[0,1],color='green',linewidth=2)
plt.xlim(xmin=0,xmax=1785)
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
nuratio1=np.zeros((len(blzdays),1))
n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.995)
    nuratio1[n]=ffrac[n]/(nulim[n])
    n+=1

nulim=np.zeros((len(blzdays),1))
nuratio2=np.zeros((len(blzdays),1))
n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.995)
    nuratio2[n]=ffrac[n]/(nulim[n]-nurate*blzdays[n])
    n+=1

nulim=np.zeros((len(blzdays),1))
nuratio3=np.zeros((len(blzdays),1))
n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.9)
    nuratio3[n]=ffrac[n]/(nulim[n]-nurate*blzdays[n])
    n+=1

nulim=np.zeros((len(blzdays),1))
nuratio4=np.zeros((len(blzdays),1))
n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.85)
    nuratio4[n]=ffrac[n]/(nulim[n]-nurate*blzdays[n])
    n+=1

nulim=np.zeros((len(blzdays),1))
nuratio5=np.zeros((len(blzdays),1))
n=0
while n<len(blzdays):
    nulim[n]=poilim(nurate*blzdays[n],0.8)
    nuratio5[n]=ffrac[n]/(nulim[n]-nurate*blzdays[n])
    n+=1



plt.figure(2)
plt.plot(blzdays,nuratio1,color='blue')
plt.plot(blzdays,nuratio2,color='red')
plt.plot(blzdays,nuratio3,color='green')
plt.plot(blzdays,nuratio4,color='yellow')
plt.plot(blzdays,nuratio5,color='black')
plt.show()



















