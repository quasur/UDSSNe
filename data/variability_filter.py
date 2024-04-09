import numpy as np
from welch_stetson import welch_stetson
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

k_x = np.loadtxt("K_band_x-axis.txt")
j_x = np.loadtxt("J_band_x-axis.txt")

file = "month_lightcurves_J_and_K.fits"
curves = fits.open(file)
fits_data = Table(curves[1].data)
columns = np.array(curves[1].columns.names[3:])

data = np.array(fits_data)

data = np.array(data.tolist())

# Loop through all data
len = data.shape[0]
count = 0
_filter = np.full(len, False)
for i in range(len):

    k_y = data[i, 3:3 + (k_x.size * 2):2]
    k_y_err = data[i, 4:4 + (k_x.size * 2):2]

    j_y = data[i, 3 + (k_x.size * 2)::2]
    j_y_err = data[i, 4 + (k_x.size * 2)::2]

    variability = welch_stetson(k_y, k_y_err, j_y, j_y_err)

    if variability >= 0.85:
        _filter[i] = True
        count += 1

data = data[_filter, :]
print(count)
#np.savetxt("wst_variable_objects.txt", data)