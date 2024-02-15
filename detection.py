#%%
import numpy as np
from scipy.stats import chisquare
#import dataset

kdat = np.loadtxt("data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("data/jdata.npy").reshape((114243,35,3))


kydat = np.loadtxt("data/kydata.npy").reshape((114243,7,3))
jydat = np.loadtxt("data/jydata.npy").reshape((114243,8,3))

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
    
plt.plot(np.log(wst))

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