import numpy as np
import matplotlib.pyplot as plt
from chi_square import chi_square
from time_series_generator import generate_time_series
import scipy.stats

steps = 20
total_runs = 10
z_runs = np.linspace(0.5, 3, steps)
semester_runs = np.linspace(0, 80, steps)
m_runs = (-22, -23, steps)

count = 0
variability = np.zeros((steps, steps))
for run in range(total_runs):
    print(run)
    for i, z in enumerate(z_runs):
        for j, ps in enumerate(semester_runs):
             #for k, m in enumerate(m_runs):
             k_x, k_y, k_y_err, j_x, j_y, j_y_err = generate_time_series(z, ps, -24)
             k_chi = chi_square(k_y, k_y_err)
             j_chi = chi_square(j_y, j_y_err)

             #k_1sig = scipy.stats.chi2.ppf(0.6825, df=k_x.size-1)
             #k_2sig = scipy.stats.chi2.ppf(0.9545, df=k_x.size - 1)
             #k_3sig = scipy.stats.chi2.ppf(0.9973, df=k_x.size - 1)
             #k_4sig = scipy.stats.chi2.ppf(0.999937, df=k_x.size - 1)
             variability[i][j] += (k_chi > 90)
             count += (k_chi > 90)
                # variability[i, j] += (k_chi > 90)
                # variability[i, j] += (k_chi > k_2sig)
                # variability[i, j] += (k_chi > k_3sig)
                # variability[i, j] += (k_chi > k_4sig)

        #print(f"run: {run+1}, i: {i}")

# np.savetxt("significance_heatmap_-23.5_expanded", variability)

# variability = np.array([np.loadtxt("significance_heatmap_-23.5_expanded"), np.loadtxt(
#     "significance_heatmap_-23_expanded"), np.loadtxt("significance_heatmap_-22.5_expanded")])
# labels = ["M=-23.5", "M=-23.0", "M=-22.5"]
#
# z_arr = np.zeros(steps)
# sig_avg_arr = np.zeros(steps)
# for j in range(3):
#     for i, z in enumerate(z_runs):
#         sig_below_z = variability[j, i:i+1, :]
#         sig_avg_arr[i] = np.mean(sig_below_z)
#         z_arr[i] = z
#
#     plt.scatter(z_arr, sig_avg_arr, label = labels[j])
# plt.xlabel("z")
# plt.ylabel("Average Detection Significance")
# plt.legend()
# plt.axhline(3, c="black", linestyle="--")
    #plt.figure()

plt.imshow(variability/total_runs, cmap="winter", extent=[0, 80, 3, 0.5], aspect='auto')
plt.xlabel("Peak Semester")
plt.ylabel("z")
plt.colorbar()
plt.show()

print(count/(steps*steps*total_runs) * 100)