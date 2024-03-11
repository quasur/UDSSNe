#%%
import numpy as np
import matplotlib.pyplot as plt
from chi_square import chi_square
from time_series_generator import generate_time_series

steps = 50
z_runs = np.linspace(1, 3, steps)
semester_runs = np.linspace(0, 80, steps)

variability1 = np.zeros((steps, steps))
variability2 = np.zeros((steps, steps))

#welch stetson test
def delta(data,error):
    mean = np.mean(data)
    return (mean - data)/error

def WST(k_y,j_y,k_y_err,j_y_err):
    #wst = np.zeros(np.size(k_y)) 
    #nudge = 0
    kdelta = delta(k_y,k_y_err)
    jdelta = delta(j_y  ,j_y_err)
    new = np.sqrt(1/(35*34))*np.sum(kdelta*jdelta)
    return new


for i, z in enumerate(z_runs):
    for j, ps in enumerate(semester_runs):
        k_x, k_y, k_y_err, j_x, j_y, j_y_err = generate_time_series(z, ps)
        k_chi = chi_square(k_y, k_y_err)
        j_chi = chi_square(j_y, j_y_err)

        welchst = WST(k_y[0:35],j_y,k_y_err[0:35],j_y_err)


        variability1[i, j] = welchst
        variability2[i,j]  = k_chi


    print(i)
#%%
plt.imshow(variability2, cmap="plasma", extent=[0, 80, 3, 1], aspect='auto')
plt.colorbar()
plt.show()

