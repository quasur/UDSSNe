import numpy as np
import time_series_generator as tsg
import background_generator as bg
from chi_square import chi_square
import scipy.stats
import matplotlib.pyplot as plt

def get_fpr(z, magnitude, magnitude_range, runs):

    # Get data from UDS that has magnitudes in range
    k_y_err_mean, j_y_err_mean = bg.get_data_in_mag_range(z, magnitude, magnitude_range)

    k_1sig = scipy.stats.chi2.ppf(0.6825, df=38 - 1)
    k_2sig = scipy.stats.chi2.ppf(0.9545, df=38 - 1)
    k_3sig = scipy.stats.chi2.ppf(0.9973, df=38 - 1)
    k_4sig = scipy.stats.chi2.ppf(0.999937, df=38 - 1)

    variance = np.zeros(runs)
    for i in range(runs):
        # Generate background
        background_k, background_k_err, background_j, background_j_err = bg.generate_background_from_data(k_y_err_mean, j_y_err_mean)
        #background_k = 1
        #background_k_err = 1

        background_k = np.zeros(38)
        background_k += np.random.normal(0, background_k_err)
        k_chi = chi_square(background_k, background_k_err)

        variance[i] += (k_chi > k_1sig)
        variance[i] += (k_chi > k_2sig)
        variance[i] += (k_chi > k_3sig)
        variance[i] += (k_chi > k_4sig)

    return np.count_nonzero(variance >= 3) / runs * 100


def get_tpr(z, magnitude, magnitude_range, runs):

    k_x = np.loadtxt("K_band_x-axis.txt")
    j_x = np.loadtxt("J_band_x-axis.txt")

    # Get data from UDS that has magnitudes in range
    k_y_err_mean, j_y_err_mean = bg.get_data_in_mag_range(z, magnitude, magnitude_range)

    # Get semester range
    semester_range = np.linspace(0, 80, runs)

    k_1sig = scipy.stats.chi2.ppf(0.6825, df=38 - 1)
    k_2sig = scipy.stats.chi2.ppf(0.9545, df=38 - 1)
    k_3sig = scipy.stats.chi2.ppf(0.9973, df=38 - 1)
    k_4sig = scipy.stats.chi2.ppf(0.999937, df=38 - 1)

    variance = np.zeros(runs)
    for i in range(runs):
        print(i)
        background_k, background_k_err, background_j, background_j_err = bg.generate_background_from_data(k_y_err_mean, j_y_err_mean)
        k_x, k_y, k_y_err, j_x, j_y, j_y_err = tsg.generate_time_series(z, semester_range[i], magnitude, k_x, background_k, background_k_err, j_x, background_j, background_j_err)
        #k_y = 1
        #k_y_err = 1
        k_chi = chi_square(k_y, k_y_err)

        variance[i] += (k_chi > k_1sig)
        variance[i] += (k_chi > k_2sig)
        variance[i] += (k_chi > k_3sig)
        variance[i] += (k_chi > k_4sig)


    return np.count_nonzero(variance >= 3) / runs * 100


redshift = np.linspace(0.5, 3, 10)
tpr_arr = np.zeros(redshift.size)
fpr_arr = np.zeros(redshift.size)
mags = np.array([-22.5, -23, -23.5])

for j, m in enumerate(mags):
    for i, z in enumerate(redshift):
        print(z)
        #tpr = get_tpr(z, -23, 0.5, 10000)
        fpr = get_fpr(z, mags[j], 0.5, 100000)

        #tpr_arr[i] = tpr
        fpr_arr[i] = fpr

    plt.scatter(redshift, fpr_arr, label = mags[j].astype(str))
#plt.plot(redshift, fpr_arr)
plt.legend()
plt.show()