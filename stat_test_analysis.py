import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare
from curve_generator import generate_time_series


steps = 50
z_runs = np.linspace(1, 3, steps)
semester_runs = np.linspace(0, 80, steps)

variability = np.zeros((steps, steps))

for i, z in enumerate(z_runs):
    for j, ps in enumerate(semester_runs):
        k_x, k_y, k_y_err, j_x, j_y, j_y_err = generate_time_series("Nugent Templates\\sn1a_flux.v1.2.dat", z, ps)
        k_chi = chisquare(k_y, k_y_err)
        j_chi = chisquare(j_y, j_y_err)

        variability[i, j] = k_chi

    print(i)

plt.imshow(variability, cmap="plasma", extent=[0, 80, 3, 1], aspect='auto')
plt.colorbar()
plt.show()


"""
#welch stetson test
def delta(data,error):
    mean = np.mean(data)
    return (mean - data)/error

wst = np.zeros(np.size(kdat[:,0,0])) 
nudge = 0
for j in range(nudge+3):
    for i in range(np.size(kdat[:,0,0])):
        kdelta = delta(kdat[i,j:j+35-nudge,1],kdat[i,j:j+35-nudge,2])
        jdelta = delta(jdat[i,0:35-nudge,1]  ,jdat[i,0:35-nudge,2])
        new = np.sqrt(1/(35*34))*np.sum(kdelta*jdelta)
        if new > wst[i] or wst[i]==0:
            wst[i]=new
"""
