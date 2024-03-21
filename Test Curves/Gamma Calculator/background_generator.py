import numpy as np
from astropy.io import fits
from astropy.table import Table
from luminosity_distance import redshift_to_lum_distance
import scipy.stats
from chi_square import chi_square
import matplotlib.pyplot as plt

def get_data_in_mag_range(z, magnitude, mag_range):

    k_x = np.loadtxt("..\\SN Curves\\K_band_x-axis.txt")

    # Import fits data
    file = "..\\..\\data\\month_lightcurves_J_and_K.fits"
    curves = fits.open(file)
    fits_data = Table(curves[1].data)

    df = fits_data.to_pandas()
    data = df.to_numpy()

    magnitude = magnitude
    redshift = z
    dist = redshift_to_lum_distance(redshift)
    apparent_mag = 5*np.log10(dist * 1e6/10) + magnitude

    # Get average at uncertainty for magnitude
    mag_filter = np.zeros(data.shape[0])
    for i in range(data.shape[0]):

        k_y_avg = np.mean(data[i, 3:3 + (k_x.size * 2):2], axis=0)

        apparent_mag_DR11 = 30-2.5*np.log10(k_y_avg)

        mag_filter[i] = (np.abs(apparent_mag_DR11-apparent_mag) < mag_range)

    k_y_err = data[:, 4:4 + (38 * 2):2]
    k_y_err_mean = np.mean(k_y_err, axis=0)

    j_y_err = data[:, 4 + (38 * 2)::2]
    j_y_err_mean = np.mean(j_y_err, axis=0)

    return k_y_err_mean, j_y_err_mean




def generate_background_from_data(k_y_err_mean, j_y_err_mean):

    k_y_err_mean = k_y_err_mean
    j_y_err_mean = j_y_err_mean
    # Get errors from UDS data

    # K band
    k_y = np.zeros(38)
    k_y += np.random.normal(0, k_y_err_mean)

    # J band
    j_y = np.zeros(35)
    j_y += np.random.normal(0, j_y_err_mean)


    return k_y, k_y_err_mean, j_y, j_y_err_mean



#
# def get_FPR(z, magnitude, magnitude_range, runs):
#     k_x = np.loadtxt("..\\SN Curves\\K_band_x-axis.txt")
#
#     data = get_data_in_mag_range(z, magnitude, magnitude_range)
#     # Get errors from UDS data
#     k_y_err = data[:, 4:4 + (k_x.size*2):2]
#     k_y_err_mean = np.mean(k_y_err, axis=0)
#
#     runs = runs
#
#     variance = np.zeros(runs)
#
#     k_1sig = scipy.stats.chi2.ppf(0.6825, df=k_x.size - 1)
#     k_2sig = scipy.stats.chi2.ppf(0.9545, df=k_x.size - 1)
#     k_3sig = scipy.stats.chi2.ppf(0.9973, df=k_x.size - 1)
#     k_4sig = scipy.stats.chi2.ppf(0.999937, df=k_x.size - 1)
#
#     for i in range(runs):
#         #print(i)
#         k_y = np.zeros(k_x.size)
#         k_y += np.random.normal(0, k_y_err_mean)
#         k_chi = chi_square(k_y, k_y_err_mean)
#
#         variance[i] += (k_chi > k_1sig)
#         variance[i] += (k_chi > k_2sig)
#         variance[i] += (k_chi > k_3sig)
#         variance[i] += (k_chi > k_4sig)
#
#     return np.count_nonzero(variance >= 3)/runs * 100
#
# get_FPR(1.5, -23, 0.5, 10000)