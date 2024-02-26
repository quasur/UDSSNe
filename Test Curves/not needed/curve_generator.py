import numpy as np
import matplotlib.pyplot as plt


def generate_template(template_path, z, peak_semester, flux):
    path = template_path
    template_xdata = np.loadtxt("template_xdata.txt")
    template_ydata = np.loadtxt("template_ydata.txt")
    template_xdata /= (30 * 24)

    # Move peak to x = 0
    template_xdata -= template_xdata[np.argmax(template_ydata)]

    # Stretch based on redshift
    template_xdata *= (1 + z)

    # Move peak to peak_semester
    template_xdata += peak_semester

    template_ydata *= ((1/template_ydata.max()) * flux)

    return template_xdata, template_ydata
