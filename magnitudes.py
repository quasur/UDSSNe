from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from luminosity_distance import redshift_to_lum_distance

ids = np.array([153802, 111791, 62253, 104987, 147703])
i = ids[4]

hdul = fits.open("C:\\Users\\jcamp\\OneDrive\\Desktop\\University\\Year 4\\Project\\Light Curves\\DR11.fits")
redshifts = np.array(hdul[1].data['z_p'])
r_ids = np.array(hdul[1].data['ID'])

K_x = np.loadtxt("data\\K_band_x-axis.txt")
J_x = np.loadtxt("data\\J_band_x-axis.txt")
hdul = fits.open("data\\month_lightcurves_J_and_K.fits")
curve_data = Table(hdul[1].data)

dr11_ids = np.array(hdul[1].data['DR11_ID'])

curve_index = np.where(dr11_ids == i)[0][0]
redshift_index = np.where(r_ids == i)[0][0]

data = np.array(curve_data)[curve_index]
data = np.array(data.tolist())

K_y = data[3:3 + (K_x.size*2):2]
K_y_err = data[4:4 + (K_x.size*2):2]

J_y = data[3 + (K_x.size*2)::2]
J_y_err = data[4 + (K_x.size*2)::2]

dist = redshift_to_lum_distance(redshifts[redshift_index])

F_k = K_y.max() - K_y.mean()
F_10_k = np.square(1e5*dist)*F_k

F_j = J_y.max() - J_y.mean()
F_10_j = np.square(1e5*dist)*F_j
print(redshifts[redshift_index])
M_k = 30 - (2.5*np.log10(F_10_k)) + 1.9
M_j = 30 - (2.5*np.log10(F_10_j)) + 0.938

print(f"K-band magnitude: {M_k}, J-band magnitude: {M_j}")