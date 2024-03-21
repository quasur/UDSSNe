import numpy as np

def chi_square(ydata, ydata_err):
    variation = ((ydata - np.mean(ydata))**2)/(ydata_err**2)
    return np.sum(variation)