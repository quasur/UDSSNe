import numpy as np

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
