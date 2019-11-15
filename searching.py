import numpy as np
import matplotlib.pyplot as plt

def icscramble(data,seed):
    #first, we creat a blank array of the same dimensions as the data
    base=np.zeros((len(data[:,0]),3))
    #then we set its columns equal to the data, randomizing arrival time
    np.random.seed(seed)
    base[:,0]=np.random.permutation(data[:,0])
    base[:,1]=data[:,1]
    base[:,2]=data[:,2]
    #now we set up a loop to fix the right ascension
    n=0 #loop counter
    siderial=0.99726958 #length of a siderial day in solar days
    while n<len(data[:,1]):
        diff=base[n,0]-data[n,0]
        base[n,1]+=diff*24/siderial
        while base[n,1]>24:
            base[n,1]-=24
        while base[n,1]<0:
            base[n,1]+=24
        n+=1
    # outside the loop, we sort based on time once again     
    base=base[np.argsort(data[:,0])] 
    return base

def distsph(lng1,lat1,lng2,lat2):
    hahaha=np.arccos(np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lng1-lng2))
    return hahaha*180/np.pi

ic40=np.loadtxt('IC40.dat',usecols=(8,1,0))
ic40=ic40[np.where(ic40[:,2]>25)]


mrk421=np.loadtxt('itms/longtermpairs.txt')
blzra=mrk421[0,0]*np.pi/180
blzdec=mrk421[0,1]*np.pi/180

es0806=np.loadtxt('itms/1es0806+524ontimes.txt')
ra0=es0806[0,0]*np.pi/180
dec0=es0806[0,1]*np.pi/180

es1218=np.loadtxt('itms/1es1218+304ontimes.txt')
ra1=es1218[0,0]*np.pi/180
dec1=es1218[0,1]*np.pi/180

c66a=np.loadtxt('itms/3c66aontimes.txt')
ra2=c66a[0,0]*np.pi/180
dec2=c66a[0,1]*np.pi/180

mrk501=np.loadtxt('itms/mrk501ontimes.txt')
ra3=mrk501[0,0]*np.pi/180
dec3=mrk501[0,1]*np.pi/180

wcomae=np.loadtxt('itms/wcomaeontimes.txt')
ra4=wcomae[0,0]*np.pi/180
dec4=wcomae[0,1]*np.pi/180

count=10000
np.random.seed(321619895)
seeds=np.random.random_integers(10**6,10**9,(count,1))
n=0
hst=np.zeros((1,1))
mhst=np.zeros((1,1))
while n<count:
    scb=icscramble(ic40,seeds[n])
    icmjd=scb[:,0]
    icra=scb[:,1]*np.pi/12
    icdec=scb[:,2]*np.pi/180
    m=0
    detected=0
    mrkdet=0
    while m<len(icmjd):
        dst=distsph(blzra,blzdec,icra[m],icdec[m])
        dst0=distsph(ra0,dec0,icra[m],icdec[m])
        dst1=distsph(ra1,dec1,icra[m],icdec[m])
        dst2=distsph(ra2,dec2,icra[m],icdec[m])
        dst3=distsph(ra3,dec3,icra[m],icdec[m])
        dst4=distsph(ra4,dec4,icra[m],icdec[m])
        
        if dst<=2.3:
            k=1
            while k<len(mrk421[:,0]):
                if mrk421[k,0]<=icmjd[m]<=mrk421[k,1]:
                    mrkdet+=1
                k+=1
        if dst0<2.3:
            if es0806[1,0]<=icmjd[m]<=es0806[1,1]:
                detected+=1
        
        elif dst1<=2.3:
            k=1
            while k<len(es1218[:,0]):
                if es1218[k,0]<=icmjd[m]<=es1218[k,1]:
                    detected+=1
                k+=1
                
        elif dst2<=2.3:
            k=1
            while k<len(c66a[:,0]):
                if c66a[k,0]<=icmjd[m]<=c66a[k,1]:
                    detected+=1
                k+=1
            
        elif dst3<=2.3:
            k=1
            while k<len(mrk501[:,0]):
                if mrk501[k,0]<=icmjd[m]<=mrk501[k,1]:
                    detected+=1
                k+=1
                
        elif dst4<=2.3:
            k=1
            while k<len(wcomae[:,0]):
                if wcomae[k,0]<=icmjd[m]<=wcomae[k,1]:
                    detected+=1
                k+=1
                
                
        m+=1
    holder1=np.zeros((1,1))
    holder2=np.zeros((1,1))
    holder1[0,0]=mrkdet
    holder2[0,0]=detected
    mhst=np.append(mhst,holder1,0)
    hst=np.append(hst,holder2,0)      
    n+=1
    if n%10==0:
        print n
    
    
    
mhst=mhst[1:len(hst),:]
hst=hst[1:len(hst),:]


plt.figure(1)
plt.title('mrk 421')
plt.hist(mhst,bins=[0,1,2,3,4,5,6,7,8])
plt.show()

plt.figure(2)
plt.title('others')
plt.hist(hst,bins=[0,1,2,3,4,5,6,7,8,9,10])
plt.show()

time1=np.zeros((1,2))
time2=np.zeros((1,2))
time1[0,0]=np.argmax(mhst)
time1[0,1]=max(mhst)
time2[0,0]=np.argmax(hst)
time2[0,1]=max(hst)




