import numpy as np
import matplotlib.pyplot as plt
import luminosity_distance


def generate_template(z, peak_semester):
    template_xdata = np.loadtxt("template_xdata.txt")
    template_ydata = np.loadtxt("template_ydata.txt")
    template_xdata /= (30 * 24)

    # Move peak to x = 0
    template_xdata -= template_xdata[np.argmax(template_ydata)]

    # Stretch based on redshift
    template_xdata *= (1 + z)

    # Move peak to peak_semester
    template_xdata += peak_semester

    template_ydata *= (1 / template_ydata.max())

    return template_xdata, template_ydata


def generate_time_series(z, peak_semester, magnitude, plot=False):
    k_x = np.loadtxt("K_band_x-axis.txt")
    j_x = np.loadtxt("J_band_x-axis.txt")

    flux_arr = np.loadtxt("background_flux_library.txt")
    background_index = np.random.randint(flux_arr.shape[0])

    k_y = flux_arr[background_index][3:3 + (k_x.size * 2):2]
    k_y_err = flux_arr[background_index][4:4 + (k_x.size * 2):2]

    j_y = flux_arr[background_index][3 + (k_x.size * 2):3 + (k_x.size * 4):2]
    j_y_err = flux_arr[background_index][4 + (k_x.size * 2)::2]

    k_y_bg = k_y
    j_y_bg = j_y

    # Get luminosity distance from redshift
    ld = luminosity_distance.redshift_to_lum_distance(z)
    ld_known = luminosity_distance.redshift_to_lum_distance(1.5)

    k_flux_diff = 10**((30-magnitude+1.9)/2.5) * (10/(ld*1e6))**2
    j_flux_diff = 10 ** ((30 - magnitude+0.938) / 2.5) * (10 / (ld * 1e6)) ** 2

    # Get flux from luminosity distance relationship
    sn_flux_k = k_flux_diff * (ld_known / ld) ** 2
    sn_flux_j = j_flux_diff * (ld_known / ld) ** 2

    # Get the x and y arrays for the template, given a peak_semester and redshift
    dummy_curve_x, dummy_curve_y = generate_template(z, peak_semester)

    k_y = k_y_bg - np.mean(k_y_bg) + np.interp(k_x, dummy_curve_x, dummy_curve_y * sn_flux_k)
    j_y = j_y_bg - np.mean(j_y_bg) + np.interp(j_x, dummy_curve_x, dummy_curve_y * sn_flux_j)

    if plot == True:
        fig = plt.figure(constrained_layout=True)
        label_ax = fig.add_subplot(111)
        ax = [fig.add_subplot(211), fig.add_subplot(212)]

        label_ax.spines['top'].set_color('none')
        label_ax.spines['bottom'].set_color('none')
        label_ax.spines['left'].set_color('none')
        label_ax.spines['right'].set_color('none')
        label_ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

        ax[0].errorbar(k_x, k_y, yerr=k_y_err, fmt='o', label="K-band")
        ax[0].errorbar(j_x, j_y, yerr=j_y_err, fmt='o', label="J-band")

        ax[1].plot(np.append(np.insert(dummy_curve_x, 0, -50), 120),
                   np.append(np.insert(dummy_curve_y * sn_flux_k, 0, 0), 0), label="K-band")
        ax[1].plot(np.append(np.insert(dummy_curve_x, 0, -50), 120),
                   np.append(np.insert(dummy_curve_y * sn_flux_j, 0, 0), 0), label="J-band")

        ax[0].tick_params(axis='x', direction='inout')

        ax[1].set_xlim(ax[0].get_xlim())
        ax[1].tick_params(axis='x', direction='inout')
        ax[1].xaxis.set_ticks_position('both')

        ax[1].set_xlabel("Month since initial observation", labelpad=15)
        ax[0].set_ylabel(r"$\mathcal{F} - \bar{\mathcal{F}}_{bg}$", labelpad=7, fontsize=15)
        ax[1].set_ylabel(r"$\mathcal{F}$", labelpad=10, fontsize=15)

        ax[0].legend()
        ax[1].legend()

        fig.subplots_adjust(hspace=0, wspace=0)
        plt.show()

    return k_x, k_y, k_y_err, j_x, j_y, j_y_err



#kx, K_y, k_y_err, jx, J_y, j_y_err = generate_time_series(2, 50, -23, plot=True)