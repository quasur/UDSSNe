import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from scipy import stats
import statsmodels.api as sm
import os


def mean_(val, freq):
    return np.average(val, weights=freq)


def var_(val, freq):
    avg = mean_(val, freq)
    dev = freq * (val - avg) ** 2
    return dev.sum() / (freq.sum() - 1)


def std_(val, freq):
    return np.sqrt(var_(val, freq))


# Import yearly data
yr_hdul = fits.open("C:\\Users\\jcamp\\OneDrive\\Desktop\\University\\Year 4\\Project\\Light Curves\\semester_lightcurves_J_and_K.fits")
yr_curves = Table(yr_hdul[1].data)
yr_curves = yr_curves.to_pandas()
yr_curves = yr_curves.to_numpy()
yr_curve_ids = np.array(yr_hdul[1].data['DR11_ID'])

# Import image data
objects = np.loadtxt("data\\xray_filtered_variables.txt")
object_ids = objects[:, 0]
print(np.where(object_ids == 22))
# Open fits file with object image coordinates
#with fits.open("C:\\Users\\jcamp\\OneDrive\\Desktop\\University\\Year 4\\Project\\Light Curves\\DR11.fits") as hdul:  # open a FITS file
#    dr11_data = hdul[1].data
hdul = fits.open("C:\\Users\\jcamp\\OneDrive\\Desktop\\University\\Year 4\\Project\\Light Curves\\DR11.fits")
dr11_ids = np.array(hdul[1].data['ID'])
dr11_x = np.array(hdul[1].data['X_IMAGE'])
dr11_y = np.array(hdul[1].data['Y_IMAGE'])

sn_detection = np.full(object_ids.size, False)
for j, oid in enumerate(object_ids):
    oid = 252582
    # Find image coordinate for each id
    array_index = np.where(dr11_ids == oid)[0][0]

    yr_index = np.where(yr_curve_ids == oid)[0][0]
    yr_curve_data = np.array(yr_curves)[yr_index]

    k_x = np.array([0, 2, 3, 4, 5, 6, 7])
    j_x = np.array([0, 1, 2, 3, 4, 5, 6, 7])

    k_y = yr_curve_data[3:3 + (k_x.size*2):2]
    k_y_err = yr_curve_data[4:4 + (k_x.size*2):2]

    j_y = yr_curve_data[3 + (k_x.size*2)::2]
    j_y_err = yr_curve_data[4 + (k_x.size*2)::2]

    y = dr11_x[array_index]
    x = dr11_y[array_index]
    frame_size = 20
    x_bounds = (int(x - frame_size), int(x + frame_size))
    y_bounds = (int(y - frame_size), int(y + frame_size))

    ax = np.array([])
    ax = np.append(ax, plt.subplot2grid((3, 8), (0, 0), rowspan=1, colspan=8))
    ax = np.append(ax, plt.subplot2grid((3, 8), (1, 0), rowspan=1, colspan=8))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 0), rowspan=1, colspan=1))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 2), rowspan=1, colspan=1))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 3), rowspan=1, colspan=1))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 4), rowspan=1, colspan=1))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 5), rowspan=1, colspan=1))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 6), rowspan=1, colspan=1))
    ax = np.append(ax, plt.subplot2grid((3, 8), (2, 7), rowspan=1, colspan=1))

    for i in range(ax.size-2):
        ax[i+2].axis('off')

    # Plot curves
    ax[0].errorbar(k_x, k_y, yerr=k_y_err, fmt='o', c="green", label="K-band")
    ax[1].errorbar(j_x, j_y, yerr=j_y_err, fmt='o', c="blue", label="J-band")
    #ax[0].text(0.01, 1.07, f"DR11 ID: {int(oid)}", ha='left', va='top', transform=ax[0].transAxes)
    ax[0].set_title(f"DR11 ID: {int(oid)}")
    ax[0].legend()
    ax[1].legend()

    # Import images
    directory = "D:\\sem_imgs"
    data = np.zeros((7, 2 * frame_size, 2 * frame_size))
    counter = 0
    for filename in os.listdir(directory):

        # Check the file is a fits file
        if filename.endswith(".fits"):
            file = os.path.join(directory, filename)

            if "06B" in filename:
                continue

            with fits.open(file) as hdul:  # open a FITS file
                image = hdul[0].data

            data[counter, :, :] = image[x_bounds[0]:x_bounds[1], y_bounds[0]:y_bounds[1]]

            counter += 1

    # Decompose into x and y distributions
    vals = np.linspace(0, 2 * frame_size, 2 * frame_size, endpoint=True)

    counts = np.zeros(data.shape[0])
    for i in range(data.shape[0]):

        ax[i + 2].imshow(data[i])
        #ax[i + 1].contour(data[i], 2, cmap="plasma")
        x_dist = np.sum(data[i], axis=0)
        y_dist = np.sum(data[i], axis=1)

        counts[i] = np.count_nonzero(data[i] > data[i].max()/1.5)

        mean = mean_(vals, x_dist)
        std = std_(vals, x_dist)
        #print(f"x mean: {mean}, x standard_deviation: {std}")

        #print(ks_x)

        # plt.figure()
        # ax1 = plt.subplot2grid((4, 4), (1, 0), colspan=3, rowspan=3)
        # ax1.pcolormesh(data[i])
        # ax1.contour(data[i], 3, cmap="plasma")
        #
        # ax2 = plt.subplot2grid((4, 4), (1, 3), rowspan=3, sharey=ax1)
        # ax2.plot(y_dist, range(len(y_dist)))
        # ax2.axis('on')
        #
        # ax3 = plt.subplot2grid((4, 4), (0, 0), colspan=3, sharex=ax1)
        # ax3.plot(range(len(x_dist)), x_dist)
        # #ax3.hist(np.random.normal(mean, std, 60000), bins=2*frame_size)
        # ax3.axis('on')
        #
        # plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()


    # counts_outside_2sigma = abs(counts-np.mean(counts)) > 2*np.std(counts)
    #
    # if np.count_nonzero(counts_outside_2sigma) == 1:
    #     sn_detection[j] = True
    # elif np.count_nonzero(counts_outside_2sigma) == 2:
    #     args = np.argwhere(counts_outside_2sigma == True)
    #     if abs(args[1]-args[0]) == 1:
    #         sn_detection[j] = True

    #print(sn_detection)

#np.savetxt("objects_to_check.txt", objects[sn_detection])


