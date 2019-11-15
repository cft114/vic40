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

def distsph(ra1,dec1,ra2,dec2):
    hahaha=np.arccos(np.sin(dec1)*np.sin(dec2)+np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2))
    return hahaha*180/np.pi


ic40=np.loadtxt('data/IC40.dat',usecols=(8,1,0,3))
ic40=ic40[np.where(ic40[:,2]>35)]

mrk421=np.loadtxt('data/longtermpairs.txt')
blzra=mrk421[0,0]*np.pi/180
blzdec=mrk421[0,1]*np.pi/180

count=10000
np.random.seed(57342)
seeds=np.random.random_integers(10**6,10**9,(count,1))
n=0
hst=np.zeros((1,1))
while n<count:
    scb=icscramble(ic40,seeds[n])
#    scb=ic40
    icmjd=scb[:,0]
    icra=scb[:,1]*np.pi/12
    icdec=scb[:,2]*np.pi/180
    m=0
    detected=0
    while m<len(icmjd):
        dst=distsph(blzra,blzdec,icra[m],icdec[m])
        if dst<=2.3:
            k=1
            while k<len(mrk421[:,0]):
                if mrk421[k,0]<=icmjd[m]<=mrk421[k,1]:
                    detected+=1
                k+=1
        m+=1
    holder=np.zeros((1,1))
    holder[0,0]=detected
    hst=np.append(hst,holder,0)       
    n+=1
    if n%10==0:
        print n
hst=hst[1:len(hst),:] 

plt.figure(1)
plt.hist(hst,bins=[0,1,2,3,4,5,6,7,8])
plt.show()

a,b=np.histogram(hst,8,(0,8))

mean1=np.average(np.arange(0,8),weights=a)

var1=(np.average((np.arange(0,8)-mean1)**2,weights=a))

print mean1,var1
