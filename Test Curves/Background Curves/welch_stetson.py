import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import math

def welch_stetson(k_y, k_y_err, j_y, j_y_err):

    n = 33

    # Trim data that doesnt line up in bands
    k_remove = np.array([9, 12, 18, 19, 32])
    _k_y = np.delete(k_y, k_remove)
    _k_y_err = np.delete(k_y_err, k_remove)

    j_remove = np.array([5, 6])
    _j_y = np.delete(j_y, j_remove)
    _j_y_err = np.delete(j_y_err, j_remove)

    k_bar = np.sum(_k_y/np.square(_k_y_err)) / np.sum(1/np.square(_k_y_err))
    j_bar = np.sum(_j_y / np.square(_j_y_err)) / np.sum(1 / np.square(_j_y_err))

    k_var = (_k_y - k_bar)/_k_y_err
    j_var = (_j_y - j_bar) / _j_y_err

    variability = np.sqrt(1/(n*(n-1)))*np.sum(k_var*j_var)

    return variability


m = np.array([-20, -21, -22, -23, -24])
chi_perc = np.array([18.975, 26.75, 36.575, 48.325, 62.925])
wst_perc = np.array([23.46, 32.33, 43.65, 56.625, 73.8])

#plt.scatter(m, chi_perc, marker="s", label=r"$\chi^2$")
#plt.scatter(m, wst_perc, marker="^")
plt.plot(m, chi_perc, linestyle="--", marker="s", label=r"$\chi^2$")
plt.plot(m, wst_perc, linestyle="--", marker="^", label=r"$I$")

plt.gca().invert_xaxis()
new_list = range(math.floor(min(m)), math.ceil(max(m))+1)
plt.xticks(new_list)

plt.xlabel("Absolute Magnitude")
plt.ylabel("Percentage of objects found")

plt.legend()
plt.show()
