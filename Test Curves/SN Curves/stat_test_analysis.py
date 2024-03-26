import numpy as np
import matplotlib.pyplot as plt
from chi_square import chi_square
from time_series_generator import generate_time_series
import scipy.stats
from welch_stetson import welch_stetson


def adjacencytest(jdata,kdata,jerr,kerr,thresh):
    length = np.size(jdata)
    jpad = np.zeros(length+2)
    kpad = jpad.copy()
    jadjnum = kpad.copy()
    kadjnum = kpad.copy()
    kpad[1:-1]=kdata.copy()
    jpad[1:-1]=jdata.copy()
    for j in range(np.size(jdata)):
        i=j+1
        jadjnum[i] = np.sum(np.abs((jpad[i]-jpad[i-1:i+1])/jerr[j])>thresh)>1
        kadjnum[i] = np.sum(np.abs((kpad[i]-kpad[i-1:i+1])/kerr[j])>thresh)>1
    return jadjnum , kadjnum



steps = 20
total_runs = 25
z_runs = np.linspace(0.5, 3, steps)
semester_runs = np.linspace(0, 80, steps)
m_runs = (-22, -23, steps)

count = 0
variability = np.zeros((steps, steps))
for run in range(total_runs):
    print(run)
    for i, z in enumerate(z_runs):
        for j, ps in enumerate(semester_runs):
            # for k, m in enumerate(m_runs):
            k_x, k_y, k_y_err, j_x, j_y, j_y_err = generate_time_series(z, ps, -24)
            #var = welch_stetson(k_y, k_y_err, j_y, j_y_err)

            jadj,kadj = adjacencytest(j_y,k_y,j_y_err,k_y_err,2) 
            if np.sum(jadj)<1  and np.sum(kadj)<1:
                count += 1#(var > 0.7)
                variability[i][j] += 1#(var > 0.7)

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
# plt.figure()

plt.imshow(variability / total_runs, cmap="winter", extent=[0, 80, 3, 0.5], aspect='auto')
plt.xlabel("Peak Semester")
plt.ylabel("z")
plt.colorbar()
plt.show()

print(count / (steps * steps * total_runs) * 100)
#%%

kx,ky,kyerr,jx,jy,jyerr = generate_time_series(0.5,10,-24)
plt.plot(kx,ky)
plt.plot(jx,jy)