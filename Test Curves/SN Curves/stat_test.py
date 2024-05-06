#%%
import numpy as np
import matplotlib.pyplot as plt
import luminosity_distance

#supply the generate time series function with the 3xobjxtime array of (x,y,err)
#it will generate a background and supernova curve based on the parameters specified (z, month, mag)
#


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



#Define tests
#welch stetson test
def wstest(jdata,kdata,jerr,kerr):
    size = len(jdata)
    kdelta = (kdata-np.sum(kdata*kerr**-2)/np.sum(kerr**-2))/kerr
    jdelta = (jdata-np.sum(jdata*jerr**-2)/np.sum(jerr**-2))/jerr
    index = (1/np.sqrt(size*(size-1)))*np.sum(kdelta*jdelta)
    r=0
    if index > 0.7:
        r = 1
    return r 



#returns peak test and balance test results
def peaktest(jdata,kdata,jerr,kerr):
    javg = np.mean(jdata)
    kavg = np.mean(kdata)
    
    #deviation from mean in num of error bars
    jpeaks = (jdata-javg)/jerr
    kpeaks = (kdata-kavg)/kerr

    #jdeviations = np.abs(jpeaks)
    #kdeviations = np.abs(kpeaks)

    #if both j and k bands have a point with >4 deviations from mean then test is positive
    if np.sum(jpeaks>4)==1 and np.sum(kpeaks>4)==1:
        r=1
    else: r=0

    #if total error is greater below mean than above then test is positive. 
    if np.sum(jpeaks)<0.0 and np.sum(kpeaks)<0.0:
        s=1
    else:s=0

    return r,s

#import data
kdat = np.loadtxt("../../data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("../../data/jdata.npy").reshape((114243,35,3))


kydat = np.loadtxt("../../data/kydata.npy").reshape((114243,7,3))
jydat = np.loadtxt("../../data/jydata.npy").reshape((114243,8,3))#np.delete(,(1),axis=1)


#import redshitfts
pzraw = np.loadtxt("../../data/pz.csv",delimiter = ",")
pzbool = (pzraw[:,1] >-80)
pz = pzraw[pzbool,:]

#import lookup table
lut = np.loadtxt("../../data/LUT.npy")

#convert years to months
kydat[:,:,0]=(kydat[:,:,0]-2005)*12
jydat[:,:,0]=(jydat[:,:,0]-2005)*12

id=95778
snjy,snky,bgjy,bgky = generate_time_series(jydat[id,:,:],kydat[id,:,:],1,25,-21)


#Define tests
def wstest(jdata,kdata,jerr,kerr):
    size = len(jdata)
    kdelta = (kdata-np.sum(kdata*kerr**-2)/np.sum(kerr**-2))/kerr
    jdelta = (jdata-np.sum(jdata*jerr**-2)/np.sum(jerr**-2))/jerr
    index = (1/np.sqrt(size*(size-1)))*np.sum(kdelta*jdelta)
    r=0
    if index > 0.7:
        r = 1
    return r 


#returns the number of stds a point is from the mean
def peaktest(jdata,kdata,jerr,kerr):
    javg = np.mean(jdata)
    kavg = np.mean(kdata)
    
    #deviation from mean in num of error bars
    jpeaks = (jdata-javg)/jerr
    kpeaks = (kdata-kavg)/kerr

    #jdeviations = np.abs(jpeaks)
    #kdeviations = np.abs(kpeaks)

    #if both j and k bands have a point with >4 deviations from mean then test is positive
    if np.sum(jpeaks>4)==1 and np.sum(kpeaks>4)==1:
        r=1
    else: r=0

    #if total error is greater below mean than above then test is positive. 
    if np.sum(jpeaks)<-1 and np.sum(kpeaks)<-1:
        s=1
    else:s=0

    return r,s

#Returns the number of points adjacent to it within 3stds will be 1 if none match becase its matching itself
def adjacencytest(jdata,kdata,jerr,kerr):
    length = np.size(jdata)


    jpad = np.ones(length+2)*np.mean(jdata)
    kpad = np.ones(length+2)*np.mean(kdata)
    jadjnum = np.ones(length)
    kadjnum = np.ones(length)
    
    kpad[1:-1]=kdata
    jpad[1:-1]=jdata
    
    for j in range(length):
        i=j+1
        #if the sum of erros of the differences between ajacent points indicates a target point is >3 error from them then
        jadjnum[j] = np.sum(((jpad[i]-jpad[i-1:i+2])/jerr[j])>3)==2
        kadjnum[j] = np.sum(((kpad[i]-kpad[i-1:i+2])/kerr[j])>3)==2

    if np.sum(jadjnum) ==1 and np.sum(kadjnum)==1:
        r= True
    else:
        r=False
    return r


#wst but removing a single point
def wstestremove(jdata,jerr,kdata,kerr):
    size = np.size(jdata)
    wst = np.zeros(size)
    for i in range(size):
        currj = np.delete(jdata,i)
        currk = np.delete(kdata,i)
        currjerr = np.delete(jerr,i)
        currkerr = np.delete(kerr,i)
        wst[i]=wstest(currj,currjerr,currk,currkerr)

    if np.sum(wst<1)==1: 
        r = 1
    else: r=0
    return r

#parameters and var setup
stop = kdat[0,-1,0]
zsize=20
tsize=20
redshifts = np.linspace(0.1,2.5,zsize)
peakmonth = np.linspace(0,stop,tsize)
repeats = 50

wstSNresultarr = np.zeros((zsize,tsize))
wstBGresultarr = np.zeros((zsize,tsize))

peakSNresultarr = np.zeros((zsize,tsize))
peakBGresultarr = np.zeros((zsize,tsize))

balSNresultarr = np.zeros((zsize,tsize))
balBGresultarr = np.zeros((zsize,tsize))

adjSNresultarr = np.zeros((zsize,tsize))
adjBGresultarr = np.zeros((zsize,tsize))

wstremSNresultarr = np.zeros((zsize,tsize))
wstremBGresultarr = np.zeros((zsize,tsize))

#test testing over a range of parameters
#currently using a 
for i in range(zsize):

    #zrange = (redshifts[1]-redshifts[0])/2
    #print(redshifts[i]-zrange,redshifts[i]+zrange)
    #zlower = pz[:,1]>redshifts[i]-zrange
    #zupper = pz[:,1]<redshifts[i]+zrange
    #pzRangeBool = zlower & zupper
    #filter redshift list to range
    #zobjlist = pz[pzRangeBool,:]
    #number of objects in the range
    #numObj = np.sum(pzRangeBool)
    
    for j in range(tsize):
        blockcount = 0
        for k in range(repeats):
            
            objid = np.random.randint(114243)

            #zobjid = pz[np.random.randint(0,numObj),0]
            #objid = lut[:,0]==zobjid
            #print(np.sum(objid))
            
            samplej = jdat[objid,:,:]
            samplek = kdat[objid,:,:]
            samplejyear = np.delete(jydat[objid,:,:],1,axis=0)
            samplekyear = kydat[objid,:,:]

            snjy,snky, bgjy,bgky = generate_time_series(samplej,samplek,redshifts[i],peakmonth[j],-21)
            snjyy,snkyy, bgjyy,bgkyy = generate_time_series(samplejyear,samplekyear,redshifts[i],peakmonth[j],-21)

            if wstest(snjyy,snkyy,samplejyear[:,2],samplekyear[:,2])==1:
                wstSNresultarr[i,j] += 1
                r,s = peaktest(snjyy,snkyy,samplejyear[:,2],samplekyear[:,2])
                balSNresultarr[i,j]    += s
                peakSNresultarr[i,j]   += r
                adjSNresultarr[i,j]    += adjacencytest(snjyy,snkyy,samplejyear[:,2],samplekyear[:,2])
                wstremSNresultarr[i,j] += wstestremove(snjyy,snkyy,samplejyear[:,2],samplekyear[:,2])
                

            if wstest(bgjyy,bgkyy,samplejyear[:,2],samplekyear[:,2])==1:
                wstBGresultarr[i,j] += 1
                r,s = peaktest(bgjyy,bgkyy,samplejyear[:,2],samplekyear[:,2])
                peakBGresultarr[i,j] += r
                balSNresultarr[i,j] += s
                adjBGresultarr[i,j] += adjacencytest(snjyy,snkyy,samplejyear[:,2],samplekyear[:,2])
                wstremBGresultarr[i,j] += wstestremove(bgjyy,bgkyy,samplejyear[:,2],samplekyear[:,2])

    print(i)

    
#%%
print("WST TPR= ",np.sum(wstSNresultarr)/(zsize*tsize*repeats))
print("WST FPR= ",np.sum(wstBGresultarr)/(zsize*tsize*repeats))
print("peak TPR= ",np.sum(peakSNresultarr)/(zsize*tsize*repeats))
print("peak FPR= ",np.sum(peakBGresultarr)/(zsize*tsize*repeats))
print("bal TPR= ",np.sum(balSNresultarr)/(zsize*tsize*repeats))
print("bal FPR= ",np.sum(balBGresultarr)/(zsize*tsize*repeats))
print("adj TPR= ",np.sum(adjSNresultarr)/(zsize*tsize*repeats))
print("adj FPR= ",np.sum(adjBGresultarr)/(zsize*tsize*repeats))
print("wstrem TPR= ",np.sum(wstremSNresultarr)/(zsize*tsize*repeats))
print("wstrem FPR= ",np.sum(wstremBGresultarr)/(zsize*tsize*repeats))