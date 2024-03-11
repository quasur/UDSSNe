#%%
import numpy as np
from scipy.stats import chisquare
#import dataset

kdat = np.loadtxt("data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("data/jdata.npy").reshape((114243,35,3))


kydat = np.loadtxt("data/kydata.npy").reshape((114243,7,3))
jydat = np.loadtxt("data/jydata.npy").reshape((114243,8,3))

kydat[:,:,0]=(kydat[:,:,0]-2005)*12
jydat[:,:,0]=(jydat[:,:,0]-2005)*12

"""
chisq=np.zeros((kdat[:,0,1].size))
for i in range(kdat[:,0,0].size):
    expval = np.mean(kdat[i,:,1])
    chisq[i] = np.sum((kdat[i,:,1]-expval)**2/expval)

"""

chisq = chisquare(kydat[:,:,1],axis=1)[0]

#%%
#welch stetson test
def delta(data,error):
    mean = np.mean(data)
    return (mean - data)/error

wst = np.zeros(np.size(kdat[:,0,0])) 
nudge = 0
for j in range(nudge+3):
    for i in range(np.size(kdat[:,0,0])):
        kdelta = delta(kdat[i,j:j+35-nudge,1],kdat[i,j:j+35-nudge,2])
        jdelta = delta(jdat[i,0:35-nudge,1]  ,jdat[i,0:35-nudge,2])
        new = np.sqrt(1/(35*34))*np.sum(kdelta*jdelta)
        if new > wst[i] or wst[i]==0:
            wst[i]=new

#%%
#Varability imbalance detection

kmean = np.expand_dims(np.mean(kdat[:,:,1],axis=1),axis=-1)
jmean = np.expand_dims(np.mean(jdat[:,:,1],axis=1),axis=-1)

kUpOrDown = (kdat[:,:,1] > kmean) #create a masking array to apply to the original data to seperate into >mean and <mean
jUpOrDown = (jdat[:,:,1] > jmean)

print(np.sum(kUpOrDown)/np.size(kdat))
#this is a really stupid way of doing it but its funciton and doesnt work so I cba to change it
kdist = np.zeros(114243)
jdist = np.zeros(114243)
for i in range(114243):
    kup   = kdat[i,kUpOrDown[i,:],1]
    jup   = jdat[i,jUpOrDown[i,:],1]
    kdown = kdat[i,~kUpOrDown[i,:],1]
    jdown = jdat[i,~jUpOrDown[i,:],1]

    kuperr   = kdat[i,kUpOrDown[i,:],2]
    juperr   = jdat[i,jUpOrDown[i,:],2]
    kdownerr = kdat[i,~kUpOrDown[i,:],2]
    jdownerr = jdat[i,~jUpOrDown[i,:],2]

    kdist[i] = kdist[i] + (np.sum((kmean[i]-kup)/kuperr) + np.sum((kmean[i]-kdown)/kdownerr))/1#kmean[i]
    jdist[i] = jdist[i] + (np.sum((jmean[i]-jup)/juperr) + np.sum((jmean[i]-jdown)/jdownerr))/1#jmean[i]

#%%
import matplotlib.pyplot as plt

plt.hist(np.log(kdist),bins=1000,alpha=0.6)
plt.hist(np.log(jdist),bins=1000,alpha=0.6)
#%%
for i in range(114243):
    if kdist[i]>10:
        plt.plot(kydat[i,:,0],kydat[i,:,1],'r-')
        plt.errorbar(kdat[i,:,0],kdat[i,:,1],yerr=kdat[i,:,2],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
        print(i,kdist[i])
        plt.show()



#%%
#peak value dector
kyinterp = np.expand_dims(np.mean(kdat[:,:,1],axis=1), axis=-1)*np.ones((114243,38))
jyinterp = np.expand_dims(np.mean(jdat[:,:,1],axis=1),axis=-1)*np.ones((114243,35))
"""
for i in range(114243):
    kyinterp[i,:,1] = np.interp(kdat[i,:,0],kydat[i,:,0],kydat[i,:,1])
    jyinterp[i,:,1] = np.interp(jdat[i,:,0],jydat[i,:,0],jydat[i,:,1])
"""
    
stdvs = 4
kpeakBool = np.any((kdat[:,:,1]-kyinterp)>kdat[:,:,2]*stdvs,axis=1)
jpeakBool = np.any((jdat[:,:,1]-jyinterp)>jdat[:,:,2]*stdvs,axis=1)

peaksBool = np.logical_or(kpeakBool,jpeakBool)

kpeakObjs = kdat[peaksBool]
jpeakObjs = jdat[peaksBool]
kypeakObjs = kyinterp[peaksBool]
jypeakObjs = jyinterp[peaksBool]

import matplotlib.pyplot as plt
for id in range(np.size(kpeakObjs[:,0,0])):
    plt.plot(kpeakObjs[id,:,0],kypeakObjs[id,:],'r-')
    plt.errorbar(kpeakObjs[id,:,0],kpeakObjs[id,:,1],yerr=kpeakObjs[id,:,2],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
    
    plt.plot(jpeakObjs[id,:,0],jypeakObjs[id,:],'r-')
    plt.errorbar(jpeakObjs[id,:,0],jpeakObjs[id,:,1],yerr=jpeakObjs[id,:,2],color="green",marker="x",lw=0,elinewidth=1,capsize=1.5)
    
    plt.show()

#%%
plt.cla()
plt.hist(np.log(wst),bins=1000,color="blue",alpha=0.7)
plt.hist(np.log(chisq),bins=1000,color="red",alpha=0.7)
plt.show()

#%%

chisq = wst
import matplotlib.pyplot as plt
plt.hist((chisq))
plt.show()
plt.plot(chisq)
plt.show()

std = np.std(chisq)
mean = np.mean(chisq)

big = np.where(np.logical_and(chisq>3*std+mean,True))#,chisq<10000000*std+mean))

kdatbig = kdat[big]
jdatbig =jdat[big]

chisqlst = chisq[big]
indx = np.arange(kdat[:,0,0].size)[big]

for i in range(kdatbig[:,0,0].size):
    plt.plot(kdatbig[i,:,0],kdatbig[i,:,1],'b.')
    plt.errorbar(kdatbig[i,:,0],kdatbig[i,:,1],kdatbig[i,:,2],color='blue',linestyle='')
    
    plt.plot(jdatbig[i,:,0],jdatbig[i,:,1],'r.')
    plt.errorbar(jdatbig[i,:,0],jdatbig[i,:,1],jdatbig[i,:,2],color='red',linestyle='')

    print((chisqlst[i]-mean)/std,indx[i])
    plt.show()
lut = np.loadtxt("data/LUT.npy")
ids = lut[big]

#%%

kzero = np.any(kdat[:,:,1]==0,axis=1)
jzero = np.any(jdat[:,:,1]==0,axis=1)

posk=lut[kzero,:]
posj=lut[jzero,:]


def setMask(array): #Funcion to remove values outside of the range specified
    array= array.T
    xmin,ymin = 34.2,-5.4
    xmax,ymax = 34.8,-4.8
    mask = np.logical_or(array[:,0]<xmin,array[:,0]>xmax)
    mask |= np.logical_or(array[:,1]<ymin,array[:,1]>ymax)
    retarray = array[~mask]
    return retarray.T

xmin,ymin = 34.2,-5.4
xmax,ymax = 34.8,-4.8

kcoremask = np.logical_or(posk[:,1]<xmin,posk[:,1]>xmax)
kcoremask |= np.logical_or(posk[:,2]<ymin,posk[:,2]>ymax)

jcoremask = np.logical_or(posj[:,1]<xmin,posj[:,1]>xmax)
jcoremask |= np.logical_or(posj[:,2]<ymin,posj[:,2]>ymax)



corek=posk[~kcoremask]
corej=posj[~jcoremask]

plt.figure(dpi=120)
plt.plot(lut[:,1],lut[:,2],'b,')
plt.plot(posk[:,1],posk[:,2],'r.')
plt.plot(posj[:,1],posj[:,2],'r.')
plt.plot(corek[:,1],corek[:,2],'g.')
plt.plot(corej[:,1],corej[:,2],'g.')
plt.show()