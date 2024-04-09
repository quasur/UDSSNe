import numpy as np
from welch_stetson import welch_stetson
import matplotlib.pyplot as plt

data = np.loadtxt("wst_variable_objects.txt")
k_x = np.loadtxt("K_band_x-axis.txt")
j_x = np.loadtxt("J_band_x-axis.txt")
xray_detected = np.loadtxt("xray_detected.txt")

data = data[xray_detected.astype(bool) == False, :]

zeros = np.zeros(data.shape[0])
for i in range(zeros.size):
    k_y = data[i, 3:3 + (k_x.size * 2):2]
    j_y = data[i, 3 + (k_x.size * 2)::2]
    zeros[i] = np.logical_or(np.any(j_y == 0), np.any(k_y == 0))

data = data[zeros.astype(bool) == False, :]
print(data.shape)
np.savetxt("xray_filtered_variables.txt", data)
# print(np.any(data[:, 0] == 117014))
#
# for i in range(data.shape[0]):
#     print(data[i, 0])
#
#     k_y = data[i, 3:3 + (k_x.size * 2):2]
#     k_y_err = data[i, 4:4 + (k_x.size * 2):2]
#
#     j_y = data[i, 3 + (k_x.size * 2)::2]
#     j_y_err = data[i, 4 + (k_x.size * 2)::2]
#
#     fig, ax = plt.subplots(2, 1)
#
#     ax[0].errorbar(k_x, k_y, yerr=k_y_err, label="K-band", color="green", fmt='o')
#     ax[1].errorbar(j_x, j_y, yerr=j_y_err, label="J-band", color="blue", fmt='o')
#
#     ax[0].legend()
#     ax[1].legend()
#
#     plt.figure()
#     plt.errorbar(k_x, k_y, yerr=k_y_err, label="K-band", color="green", fmt='o')
#     plt.errorbar(j_x, j_y, yerr=j_y_err, label="J-band", color="blue", fmt='o')
#
#     plt.show()

#print(np.where(data[:, 0] == 147703))
#sn_arr = np.array([153802, 202060, 111791, 147703, 62253, 104987, 117014])
#for j in range(data.shape[0]):
    #index = np.where(data[:, 0] == sn_arr[j])[0][0]
    #print(sn_arr[j])
    #j = 1
    #j = 236

#    k_y_err = data[j, 4:4 + (k_x.size * 2):2]


#j_y_err = data[j, 4 + (k_x.size * 2)::2]

    # k_remove = np.array([9, 12, 18, 19, 32])
    # _k_y = np.delete(k_y, k_remove)
    # _k_y_err = np.delete(k_y_err, k_remove)
    #
    # j_remove = np.array([5, 6])
    # _j_y = np.delete(j_y, j_remove)
    # _j_y_err = np.delete(j_y_err, j_remove)

# x = np.linspace(0, 32, 32, endpoint=True)
# var_arr = np.zeros(32)
# for i in range(32):
#     var_arr[i] = welch_stetson(_k_y[:-(i+1)], _k_y_err[:-(i+1)], _j_y[:-(i+1)], _j_y_err[:-(i+1)])
#
# # Calculate gradient of var_arr at each point
# grad_arr = np.zeros(31)
# for i in range(31):
#     grad_arr[i] = np.abs((var_arr[i+1]-var_arr[i]) / (x[i+1] - x[i]))
#
# std = np.std(grad_arr)
# mean = np.mean(grad_arr)
# outliers = np.argwhere(grad_arr > mean + (2*std))[:, 0]
# #other_points = np.argwhere(grad_arr < mean + (1*std))[:, 0]
#
# #spike = np.logical_and(np.all(np.where(np.diff(outliers, axis=0) == 1, True, False)), outliers.size + other_points.size == 31)
#
# spike = np.all(np.where(np.diff(outliers, axis=0) == 1, True, False))
#
# if spike == True:
#     count += 1
#     print(data[index, 0])
#     # fig, ax = plt.subplots(3, 1)
#     # ax[0].plot(x, np.flip(var_arr))
#     # ax[1].plot(np.linspace(0, 31, 31, endpoint=True), np.flip(grad_arr))
#     # ax[1].axhline(mean + (0.5*std), c="black", linestyle="--")
#     # ax[2].errorbar(k_x, k_y, yerr=k_y_err, label="K-band", color="green", fmt='o')
#     # ax[2].errorbar(j_x, j_y, yerr=j_y_err, label="J-band", color="blue", fmt='o')
#     # #print(data[j, 0])
#     # fig_ = plt.figure()
#     plt.errorbar(k_x, k_y, yerr=k_y_err, label="K-band", color="green", fmt='o')
#     plt.errorbar(j_x, j_y, yerr=j_y_err, label="J-band", color="blue", fmt='o')
#
#     plt.show()