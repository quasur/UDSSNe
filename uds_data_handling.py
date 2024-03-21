import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

K_x = np.loadtxt("data\\K_band_x-axis.txt")
J_x = np.loadtxt("data\\J_band_x-axis.txt")

file = "data\\month_lightcurves_J_and_K.fits"
curves = fits.open(file)
fits_data = Table(curves[1].data)
columns = np.array(curves[1].columns.names[3:])

object_index = 95788
#object_index = 22
data = np.array(fits_data)[object_index]

data = np.array(data.tolist())
print(fits_data)
K_y = data[3:3 + (K_x.size*2):2]
K_y_err = data[4:4 + (K_x.size*2):2]

J_y = data[3 + (K_x.size*2)::2]
J_y_err = data[4 + (K_x.size*2)::2]

mean_flux = np.mean(K_y)
chi_square = np.sum(((K_y-mean_flux)**2)/(K_y_err**2))
print(f"Chi Squared: {chi_square}")

# Plot lines between adjacent points
non_adjacent_points = np.array([], dtype=int)
for i in range(K_x.size-1):
    if K_x[i+1] - 1 != K_x[i]:
        non_adjacent_points = np.append(non_adjacent_points, i+1)

K_x = np.insert(K_x, non_adjacent_points, np.nan)
K_y = np.insert(K_y, non_adjacent_points, np.nan)
K_y_err = np.insert(K_y_err, non_adjacent_points, np.nan)


# Plot lines between adjacent points
non_adjacent_points = np.array([], dtype=int)
for i in range(J_x.size-1):
    if J_x[i+1] - 1 != J_x[i]:
        non_adjacent_points = np.append(non_adjacent_points, i+1)

J_x = np.insert(J_x, non_adjacent_points, np.nan)
J_y = np.insert(J_y, non_adjacent_points, np.nan)
J_y_err = np.insert(J_y_err, non_adjacent_points, np.nan)

plt.scatter(K_x, K_y, color="green")
plt.errorbar(K_x, K_y, yerr=K_y_err, label="K-band", color="green")
plt.scatter(J_x, J_y, color="blue")
plt.errorbar(J_x, J_y, yerr=J_y_err, label="J-band", color="blue")
plt.legend()
plt.show()

