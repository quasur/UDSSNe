import numpy as np
import matplotlib.pyplot as plt
from chi_square import chi_square
from time_series_generator import generate_time_series

steps = 50
z_runs = np.linspace(1, 3, steps)
semester_runs = np.linspace(0, 80, steps)

variability = np.zeros((steps, steps))

for i, z in enumerate(z_runs):
    for j, ps in enumerate(semester_runs):
        k_x, k_y, k_y_err, j_x, j_y, j_y_err = generate_time_series("Nugent Templates\\sn1a_flux.v1.2.dat", z, ps)
        k_chi = chi_square(k_y, k_y_err)
        j_chi = chi_square(j_y, j_y_err)

        variability[i, j] = k_chi

    print(i)

plt.imshow(variability, cmap="plasma", extent=[0, 80, 3, 1], aspect='auto')
plt.colorbar()
plt.show()