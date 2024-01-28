#%%
from astropy.io import fits
from astropy.table import table
import numpy as np

f = fits.open("data/month_lightcurves_J_and_K.fits")
fitsdata = np.array(f[1].data)
messdata = np.array(fitsdata.tolist())

kdata = np.zeros((114243,38,3))
jdata = np.zeros((114243,35,3))

#format k data
for j in range(38):
    i=j
    kdata[:,j,1] = messdata[:,i*2+3]
    kdata[:,j,2] = messdata[:,i*2+4]

#format j data
for j in range(35):
    i=j+38
    jdata[:,j,1] = messdata[:,i*2+3]
    jdata[:,j,2] = messdata[:,i*2+3]

#import observation time
Jx = np.loadtxt("data/J_band_x-axis.txt")
Kx = np.loadtxt("data/K_band_x-axis.txt")

#append to main array
jdata[:,:,0]=Jx
kdata[:,:,0]=Kx

#save data as 2d array
np.savetxt("data/jdata.npy",jdata.reshape(jdata.shape[0],-1))
np.savetxt("data/kdata.npy",kdata.reshape(jdata.shape[0],-1))

#Use these lines to import jdata and kdata in other scripts:
#np.loadtxt("data/kdata.np").reshape((114243,38,3))
#np.loadtxt("data/jdata.np").reshape((114243,35,3))