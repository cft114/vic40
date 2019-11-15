import numpy as np

def distsph(lng1,lat1,lng2,lat2):
    hahaha=np.arccos(np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lng1-lng2))
    return hahaha*180/np.pi

ic40=np.loadtxt('data/IC40.dat',usecols=(8,1,0))

mrk421=np.loadtxt('itms/mrk421ontimes.txt')
blzra=mrk421[0,0]*np.pi/180
blzdec=mrk421[0,1]*np.pi/180

es0806=np.loadtxt('itms/1es1218+304ontimes.txt')
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



icmjd=ic40[:,0]
icra=ic40[:,1]*np.pi/12
icdec=ic40[:,2]*np.pi/180
m=0
detected=0
detloc=np.zeros((1,3))
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
            detloc[detected,:]=ic40[m,:]
            detected+=1
            
    
    elif dst1<=2.3:
        k=1
        while k<len(es1218[:,0]):
            if es1218[k,0]<=icmjd[m]<=es1218[k,1]:
                detloc[detected,:]=ic40[m,:]
                detected+=1
            k+=1
            
    elif dst2<=2.3:
        k=1
        while k<len(c66a[:,0]):
            if c66a[k,0]<=icmjd[m]<=c66a[k,1]:
                detloc[detected,:]=ic40[m,:]
                detected+=1
            k+=1
        
    elif dst3<=2.3:
        k=1
        while k<len(mrk501[:,0]):
            if mrk501[k,0]<=icmjd[m]<=mrk501[k,1]:
                detloc[detected,:]=ic40[m,:]
                detected+=1
            k+=1
            
    elif dst4<=2.3:
        k=1
        while k<len(wcomae[:,0]):
            if wcomae[k,0]<=icmjd[m]<=wcomae[k,1]:
                detloc[detected,:]=ic40[m,:]
                detected+=1
            k+=1
            
            
    m+=1

   
    







