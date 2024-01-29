#%%
import numpy as np
import matplotlib.pyplot as plt

#import uds data
kdat = np.loadtxt("data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("data/jdata.npy").reshape((114243,35,3))
#%%
#import sample lightcurves
sampleimp = np.loadtxt("Test Light Curve Library\\sn1bc_flux.v1.1.dat")

scale = 1e9*35000

filter = np.where(sampleimp[:,0]==7)
sample = sampleimp[filter,1:]
sample[0,:,0]=sample[0,:,0]/(30*24)
sample[0,:,1]=sample[0,:,1]*scale
plt.plot(sample[0,:,0],sample[0,:,1])
plt.show()

#%%
num=0

kinterp = np.interp(kdat[num,:,0],sample[0,:,0],sample[0,:,1])
plt.plot(kdat[num,:,0],kinterp+kdat[num,:,1],'b-')
plt.plot(kdat[num,:,0],kdat[num,:,1],'r-')


