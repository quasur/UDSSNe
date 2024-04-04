import numpy as np
import matplotlib.pyplot as plt
import luminosity_distance


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

def generate_background_from_data(jy,ky,jyerr,kyerr):

    #generate empty 
    jy = jy*0 +  np.mean(jy)
    ky = ky*0 +  np.mean(ky)
    #generate a new point based on the error of current point
    for i,val in enumerate(jyerr):
        jy[i] += np.random.normal(0, val)
    for i,val in enumerate(kyerr):
        ky[i] += np.random.normal(0, val)
    return jy, ky


def generate_time_series(jdata,kdata,z, peak_semester, magnitude):
    
    bgjy, bgky = generate_background_from_data(jdata[:,1],kdata[:,1],jdata[:,2],kdata[:,2])

    # Get luminosity distance from redshift
    ld = luminosity_distance.redshift_to_lum_distance(z)
    ld_known = 11089.08334805264 #luminosity_distance.redshift_to_lum_distance(1.5) #but as a number to save computation

    k_flux_diff = 10 ** ((30 - magnitude+ 1.9 ) / 2.5) * (10 / (ld * 1e6)) ** 2
    j_flux_diff = 10 ** ((30 - magnitude+0.938) / 2.5) * (10 / (ld * 1e6)) ** 2

    # Get flux from luminosity distance relationship
    sn_flux_k = k_flux_diff * (ld_known / ld) ** 2
    sn_flux_j = j_flux_diff * (ld_known / ld) ** 2


    # Get the x and y arrays for the template, given a peak_semester and redshift
    dummy_curve_x, dummy_curve_y = generate_template(z, peak_semester)
    
    snjy = bgjy + np.interp(jdata[:,0], dummy_curve_x, dummy_curve_y * sn_flux_j)
    snky = bgky + np.interp(kdata[:,0], dummy_curve_x, dummy_curve_y * sn_flux_k)
    
    return snjy,snky,bgjy,bgky

kdat = np.loadtxt("../../data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("../../data/jdata.npy").reshape((114243,35,3))


kydat = np.loadtxt("../../data/kydata.npy").reshape((114243,7,3))
jydat = np.loadtxt("../../data/jydata.npy").reshape((114243,8,3))#np.delete(,(1),axis=1)

#convert years to months
kydat[:,:,0]=(kydat[:,:,0]-2005)*12
jydat[:,:,0]=(jydat[:,:,0]-2005)*12

id=95778
snjy,snky,bgjy,bgky = generate_time_series(jydat[id,:,:],kydat[id,:,:],1,25,-21)

#parameters and var setup
stop = kdat[0,-1,0]
zsize=20
tsize=20
redshifts = np.linspace(0.2,4,zsize)
peakmonth = np.linspace(0,stop,tsize)
repeats = 10
SNresultarr = np.zeros((zsize,tsize))
BGresultarr = np.ones((zsize,tsize))*10

def wstest(jdata,kdata,jerr,kerr):
    size = len(jdata)
    kdelta = (kdata-np.sum(kdata*kerr**-2)/np.sum(kerr**-2))/kerr
    jdelta = (jdata-np.sum(jdata*jerr**-2)/np.sum(jerr**-2))/jerr
    index = (1/np.sqrt(size*(size-1)))*np.sum(kdelta*jdelta)
    r=0
    if index > 0.7:
        r = 1
    return r 


for i in range(zsize):
    for j in range(tsize):
        blockcount = 0
        for k in range(repeats):
            
            objid = np.random.randint(0,len(kdat[:,0,0]))
            
            samplej = np.delete(jydat[objid,:,:],(1),axis=0)
            samplek = kydat[objid,:,:]

            snjy,snky, bgjy,bgky = generate_time_series(samplej,samplek,redshifts[i],peakmonth[j],-21)
            SNresultarr[i,j] += wstest(snjy,snky,samplej[:,2],samplek[:,2])
            BGresultarr[i,j] += wstest(bgjy,bgky,samplej[:,2],samplek[:,2])
    print(i)

    
#%%
plt.imshow(SNresultarr)
plt.show()
plt.imshow(BGresultarr)
plt.show()