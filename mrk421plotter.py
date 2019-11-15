import numpy as np
import matplotlib.pyplot as plt

def icscramble(data):
    #first, we creat a blank array of the same dimensions as the data
    base=np.zeros((len(data[:,0]),3))
    #then we set its columns equal to the data, randomizing right ascension
    base[:,0]=data[:,0]
    base[:,1]=np.random.permutation(data[:,1])
    base[:,2]=data[:,2]
    #now we set up a loop to fix the timestamps 
    n=0 #loop counter
    siderial=23.9344696 #length of a siderial day in hours
    while n<len(data[:,0]):
        base[n,0]+=(data[n,1]-base[n,1])/siderial
        n+=1
    # outside the loop, we sort based on time once again     
    base=base[np.argsort(data[:,0])] 
    return base

def distsph(ra1,dec1,ra2,dec2,ctl):
    hahaha=np.sin(dec1)*np.sin(dec2)+np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
    if ctl=='deg':
        out=np.arccos(hahaha)*180/np.pi
    elif ctl=='rad':
        out=np.arccos(hahaha)
    
    return out
    
def spang(ra1,dec1,ra2,dec2,cont):
    saa=(np.cos(np.pi/2-dec2)-np.cos(distsph(ra1,dec1,ra2,dec2,'rad'))*np.cos(np.pi/2-dec1))/(np.sin(distsph(ra1,dec1,ra2,dec2,'rad'))*np.sin(np.pi/2-dec1))
    resita=np.arccos(saa)
    if ra1<ra2:
        resita=-resita
   
    if cont=='rad':
        resita=resita-np.pi/2
    elif cont=='deg':
        resita=resita*180/np.pi-90
    
    return resita

ic40=np.loadtxt('IC40.dat',usecols=(8,1,0))
mrk421=np.loadtxt('mrk421ontimes.txt')
blzdec=mrk421[0,1]*np.pi/180
blzra=mrk421[0,0]*np.pi/180


scb=icscramble(ic40)
icdec=scb[:,2]*np.pi/180
icra=scb[:,1]*np.pi/12
icmjd=scb[:,0]

n=0
count1=count2=0
act=0
allz=np.zeros((1,2))
found=np.zeros((1,2))
while n<len(icmjd):
    dst=distsph(blzra,blzdec,icra[n],icdec[n],'deg')
    ang=spang(blzra,blzdec,icra[n],icdec[n],'rad')
    if dst<8:
        allz[act,0]=dst*np.cos(ang)
        allz[act,1]=dst*np.sin(ang)
        act+=1
        allz=np.append(allz,np.zeros((1,2)),0)
    if dst<=2.3:
        count1+=1
        k=1
        while k<len(mrk421[:,0]):
            if mrk421[k,0]<=icmjd[n]<=mrk421[k,1]:
                found[count2,0]=dst*np.cos(ang)
                found[count2,1]=dst*np.sin(ang)
                found=np.append(found,np.zeros((1,2)),0)
                count2+=1
            k+=1
    n+=1
if len(found[:,0])>1:
    found=found[0:len(found[:,0])-1,:]
else:
    found=found[1:len(found[:,0]),:]
    
allz=allz[0:len(allz[:,0])-1,:]


fig=plt.figure(1)
plt.plot(allz[:,0],allz[:,1],linestyle='none',marker='o',markersize=5)
plt.plot(found[:,0],found[:,1],marker='o',markersize=7,linestyle='none',color='red')
circle1=plt.Circle((0,0),2.3,alpha=.3)
fig.gca().add_artist(circle1)

plt.xlim(xmin=-4,xmax=4)
plt.ylim(ymin=-4,ymax=4)

plt.show()













