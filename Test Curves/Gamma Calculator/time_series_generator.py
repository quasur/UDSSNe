import numpy as np
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


def generate_time_series(z, peak_semester, magnitude, k_x, background_k, background_k_err, j_x, background_j, background_j_err):

    k_y = background_k
    k_y_err = background_k_err

    j_y = background_j
    j_y_err = background_j_err

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

    return k_x, k_y, k_y_err, j_x, j_y, j_y_err

print(np.array([1])[0].astype(str))