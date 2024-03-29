import numpy as np
from astropy.io import fits
from astropy.table import Table
from luminosity_distance import redshift_to_lum_distance
import scipy.stats
from chi_square import chi_square
import matplotlib.pyplot as plt
from welch_stetson import welch_stetson


kydat = np.loadtxt("..\\..\\data\\kydata.npy").reshape((114243,7,3))
jydat = np.delete(np.loadtxt("..\\..\\data\jydata.npy").reshape((114243,8,3)),(1),axis=1)


def get_data_in_mag_range(z, magnitude, mag_range):

    k_x = np.loadtxt("..\\SN Curves\\K_band_x-axis.txt")

    # Import fits data
    file = "..\\..\\data\\semester_lightcurves_J_and_K.fits"
    curves = fits.open(file)
    fits_data = Table(curves[1].data)

    df = fits_data.to_pandas()
    data = df.to_numpy()

    magnitude = magnitude
    redshift = z
    dist = redshift_to_lum_distance(redshift)
    apparent_mag = 5*np.log10(dist * 1e6/10) + magnitude

    # Get average at uncertainty for magnitude
    mag_filter = np.zeros(data.shape[0])
    for i in range(data.shape[0]):

        k_y_avg = np.mean(data[i, 3:3 + (k_x.size * 2):2], axis=0)

        apparent_mag_DR11 = 30-2.5*np.log10(k_y_avg)

        mag_filter[i] = (np.abs(apparent_mag_DR11-apparent_mag) < mag_range)

    return data[mag_filter.astype(bool), :]


def generate_background_from_data(data):
    data = data
    # Get errors from UDS data

    # K band
    k_y_err = data[:, 4:4 + (7 * 2):2]
    k_y_err_mean = np.mean(k_y_err, axis=0)

    k_y = np.zeros(38)
    k_y += np.random.normal(0, k_y_err_mean)

    # J band
    j_y_err = data[4 + (7 * 2)::2]
    j_y_err_mean = np.mean(j_y_err, axis=0)

    j_y = np.zeros(38)
    j_y += np.random.normal(0, j_y_err_mean)

    return k_y, j_y


def get_FPR(z, magnitude, magnitude_range, runs, threshold):
    k_x = np.loadtxt("..\\SN Curves\\K_band_x-axis.txt")
    j_x = np.loadtxt("..\\SN Curves\\J_band_x-axis.txt")

    data = get_data_in_mag_range(z, magnitude, magnitude_range)
    # Get errors from UDS data
    k_y_err = data[:, 4:4 + (k_x.size*2):2]
    k_y_err_mean = np.mean(k_y_err, axis=0)

    j_y_err = data[:, 4 + (k_x.size*2)::2]
    j_y_err_mean = np.mean(j_y_err, axis=0)

    runs = runs

    #
    # k_1sig = scipy.stats.chi2.ppf(0.6825, df=k_x.size - 1)
    # k_2sig = scipy.stats.chi2.ppf(0.9545, df=k_x.size - 1)
    # k_3sig = scipy.stats.chi2.ppf(0.9973, df=k_x.size - 1)
    # k_4sig = scipy.stats.chi2.ppf(0.999937, df=k_x.size - 1)
    count = 0
    for i in range(runs):
        #print(i)
        k_y = np.zeros(k_x.size)
        k_y += np.random.normal(0, k_y_err_mean)

        j_y = np.zeros(j_x.size)

        j_y += np.random.normal(0, j_y_err_mean)

        k_chi = welch_stetson(k_y, k_y_err_mean, j_y, j_y_err_mean)
        count += (k_chi > threshold)
        #variance[i] += (k_chi > 50)
        #variance[i] += (k_chi > k_2sig)
        #variance[i] += (k_chi > k_3sig)
        #variance[i] += (k_chi > k_4sig)

    return count

def generate_template(z, peak_semester):
    template_xdata = np.loadtxt("template_xdata.txt")
    template_ydata = np.loadtxt("template_ydata.txt")
    template_xdata /= (30 * 24)

    # Move peak to x = 0
    template_xdata -= template_xdata[np.argmax(template_ydata)]

    # Stretch based on redshift
    template_xdata *= (1 + z)

    # Move peak to peak_semester
    template_xdata += peak_semester

    template_ydata *= (1 / template_ydata.max())

    return template_xdata, template_ydata

templatex = np.loadtxt("template_xdata.txt")
templatey = np.loadtxt("template_ydata.txt")

plt.plot(templatex,templatey)

kydat = np.loadtxt("..\\..\\data\\kydata.npy").reshape((114243,7,3))


def gen_year_SN(z, magnitude, magnitude_range, runs, threshold):
    k_x = np.array([0,2,3,4,5,6,7])*12#np.loadtxt("..\\SN Curves\\K_band_x-axis.txt")
    j_x = np.array([0,2,3,4,5,6,7])*12#np.loadtxt("..\\SN Curves\\J_band_x-axis.txt")
    """
    flux_arr = np.loadtxt("background_flux_library.txt")
    background_index = np.random.randint(flux_arr.shape[0])

    k_y = flux_arr[background_index][3:3 + (k_x.size * 2):2]
    k_y_err = flux_arr[background_index][4:4 + (k_x.size * 2):2]

    j_y = flux_arr[background_index][3 + (k_x.size * 2):3 + (k_x.size * 4):2]
    j_y_err = flux_arr[background_index][4 + (k_x.size * 2)::2]

    k_y_bg = k_y
    j_y_bg = j_y
    """


    ld = luminosity_distance.redshift_to_lum_distance(z)
    ld_known = luminosity_distance.redshift_to_lum_distance(1.5)

    k_flux_diff = 10**((30-magnitude+1.9)/2.5) * (10/(ld*1e6))**2
    j_flux_diff = 10 ** ((30 - magnitude+0.938) / 2.5) * (10 / (ld * 1e6)) ** 2

    # Get flux from luminosity distance relationship
    sn_flux_k = k_flux_diff * (ld_known / ld) ** 2
    sn_flux_j = j_flux_diff * (ld_known / ld) ** 2

    dummy_curve_x, dummy_curve_y = generate_template(z, peak_semester)

    k_y = k_y_bg - np.mean(k_y_bg) + np.interp(k_x, dummy_curve_x, dummy_curve_y * sn_flux_k)
    j_y = j_y_bg - np.mean(j_y_bg) + np.interp(j_x, dummy_curve_x, dummy_curve_y * sn_flux_j)

    k_y = k_y_bg - np.mean(k_y_bg) + np.interp(k_x, dummy_curve_x, dummy_curve_y * sn_flux_k)
    j_y = j_y_bg - np.mean(j_y_bg) + np.interp(j_x, dummy_curve_x, dummy_curve_y * sn_flux_j)


    data = get_data_in_mag_range(z, magnitude, magnitude_range)
    # Get errors from UDS data
    k_y_err = data[:, 4:4 + (k_x.size*2):2]
    k_y_err_mean = np.mean(k_y_err, axis=0)

    j_y_err = data[:, 4 + (k_x.size*2)::2]
    j_y_err_mean = np.mean(j_y_err, axis=0)



    return count



thresholds = np.linspace(0.6, 0.8, 10)
expected_arr = np.zeros(10)
#
for i, t in enumerate(thresholds):
    fpr = get_FPR(1.5, -22.5, 0.5, 100000, t)
    ratio = (fpr/1000000)
    expected_arr[i] = 114243*ratio
    print(f"count: {fpr}, percentage: {ratio*100}%, expected number from UDS: {114243*ratio}")

plt.plot(thresholds, expected_arr)
plt.xlabel("Threshold")
plt.ylabel("Expected False Positives")
plt.axhline(1, c="black", linestyle='--')
plt.show()

#%%

#import data we use
#kdat = np.loadtxt("..\\..\\data\\kdata.npy").reshape((114243,38,3))
#jdat = np.loadtxt("..\\..\\data\\jdata.npy").reshape((114243,35,3))

kydat = np.loadtxt("..\\..\\data\\kydata.npy").reshape((114243,7,3))
jydat = np.delete(np.loadtxt("..\\..\\data\jydata.npy").reshape((114243,8,3)),(1),axis=1)

lut = np.loadtxt("..\\..\\data\\LUT.npy")
AGNlut = np.loadtxt("..\\..\\data\\AGNLUT.npy")

#find array indexes of agn
arrAGN = AGNlut.copy()
for i in range(np.size(arrAGN)):
    arrAGN[i], = (np.where(lut[:,0].astype(int)==AGNlut[i])[0])

arrAGN = arrAGN.astype(np.int64)
#returns the number of stds a point is from the mean
def peaktest(jdata,kdata,jerr,kerr):
    javg = np.mean(jdata)
    kavg = np.mean(kdata)
    
    jpeaks = (jdata-javg)/jerr
    kpeaks = (kdata-kavg)/kerr

    jdeviations = np.abs(jpeaks)
    kdeviations = np.abs(kpeaks)

    return jpeaks,kpeaks,jdeviations,kdeviations
i=0
a,b,c,d=peaktest(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])

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
        jadjnum[i] = np.sum(((jpad[i]-jpad[i-1:i+1])/jerr[j])<thresh)==1
        kadjnum[i] = np.sum(((kpad[i]-kpad[i-1:i+1])/kerr[j])<thresh)==1
    return jadjnum , kadjnum



def get_AGNFPR(thresh):
    count = 0
    for i in arrAGN:
        jadj,kadj = adjacencytest(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2],thresh)
        if np.sum(jadj)==1  and np.sum(kadj)==1:
            count+=1
    return count

fp = np.zeros((100))
threshList = np.linspace(0,3,100)

for i in range(100):
    fp[i]=get_AGNFPR(threshList[i])

plt.plot(threshList,fp)
plt.plot([0,3],[7,7])
plt.show()

