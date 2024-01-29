#%%
import numpy as np
import matplotlib.pyplot as plt

#import uds data
kdat = np.loadtxt("data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("data/jdata.npy").reshape((114243,35,3))
#%%

#scale factor for y axis of sample
scale = 1e9*3500

#x offset for SN to appear in UDS data
offst = 20

#import sample lightcurves
sampleimp = np.loadtxt("Test Light Curve Library\\sn1bc_flux.v1.1.dat")
#seperate sample curves in dataset, can iterate over them later
filter = np.where(sampleimp[:,0]==1)

sample = sampleimp[filter,1:]
#scale data appropriately for UDS data
sample[0,:,0]=sample[0,:,0]/(30*24)+offst #convert to months and apply offset
sample[0,:,1]=sample[0,:,1]*scale

plt.plot(sample[0,:,0],sample[0,:,1])
plt.show()

num=0
#interpolate curve
kinterp = np.interp(kdat[num,:,0],sample[0,:,0],sample[0,:,1])

plt.plot(kdat[num,:,0],kinterp+kdat[num,:,1],'b-')
plt.plot(kdat[num,:,0],kdat[num,:,1],'r-')
#plt.plot(sample[0,:,0],sample[0,:,1],'g-')
plt.show()

