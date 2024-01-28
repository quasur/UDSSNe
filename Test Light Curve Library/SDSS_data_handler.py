import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

directory = "SDSS SN Curves\\data"
filename = "SN722.snphot-5.04.sum"

f = os.path.join(directory, filename)

sn_name = filename.split(".")[0]

data = np.loadtxt(f, dtype=str)
number_of_filters = 5

init_run = data[0, 13].astype(int)

# Get number of entries for each filter
filter_rows = np.zeros(number_of_filters).astype(int)
pointer = np.zeros(number_of_filters).astype(int)
for i in range(data.shape[0]):
    filter_rows[int(data[i, 2])] += 1

# Generate x and y data arrays
x_data = np.array([
    np.zeros(filter_rows[0]),
    np.zeros(filter_rows[1]),
    np.zeros(filter_rows[2]),
    np.zeros(filter_rows[3]),
    np.zeros(filter_rows[4]),
], dtype=object)

y_data = np.array([
    np.zeros(filter_rows[0]),
    np.zeros(filter_rows[1]),
    np.zeros(filter_rows[2]),
    np.zeros(filter_rows[3]),
    np.zeros(filter_rows[4]),
], dtype=object)

y_err = np.array([
    np.zeros(filter_rows[0]),
    np.zeros(filter_rows[1]),
    np.zeros(filter_rows[2]),
    np.zeros(filter_rows[3]),
    np.zeros(filter_rows[4]),
], dtype=object)

for i in range(data.shape[0]):
    _filter = (data[i, 2]).astype(int)

    x_data[_filter][pointer[_filter]] = data[i, 13].astype(int) - init_run
    y_data[_filter][pointer[_filter]] = data[i, 7].astype(float)
    y_err[_filter][pointer[_filter]] = data[i, 8].astype(float)

    pointer[_filter] += 1

# -------------------DRAWING LIGHT CURVES--------------------------

# filters = ['U', 'G', 'R', 'I', 'Z']
# markers = ['X', 'v', 's', 'o', 'D']
# colors = plt.cm.jet(np.linspace(0,1,number_of_filters))
# colors = np.roll(colors, 3, axis=0)
# plt.figure()
# for i in range(number_of_filters):
#    alpha = 0.6
#    if i == 3:
#        alpha = 1
#    plt.errorbar(x_data[i], y_data[i], yerr=y_err[i], label=f'{filters[i]} Filter', alpha=alpha, marker=markers[i], markerfacecolor='None', ls='None', color=colors[i])
# plt.ylabel("Flux")
# plt.xlabel("SDSS Run From Initial")
# plt.legend()

# plt.savefig(f"real_SN_curves\\curves\\{sn_name}", dpi=200)
# plt.close()