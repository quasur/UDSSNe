#%%
import numpy as np
import matplotlib.pyplot as plt
#import dataset

kdat = np.loadtxt("data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("data/jdata.npy").reshape((114243,35,3))


kydat = np.loadtxt("data/kydata.npy").reshape((114243,7,3))
jydat = np.delete(np.loadtxt("data/jydata.npy").reshape((114243,8,3)),(1),axis=1)

lut = np.loadtxt("data/LUT.npy")
AGNlut = np.loadtxt("data/AGNLUT.npy")

#find found SN array indexes and DR11 ids
foundSN = np.array([147703,89852,62253,66612,24605,219616,148305,214136,210486,154039,224386,217727,117014,252582,248061,153802,111791,104987])
arrSN=foundSN.copy()
for i in range(np.size(foundSN)):
    arrSN[i], = np.where(lut[:,0].astype(int)==foundSN[i])[0]

#convert agn dr11s to array indexes
arrAGN = AGNlut.copy()
for i in range(np.size(arrAGN)):
    arrAGN[i], = np.where(lut[:,0].astype(int)==AGNlut[i])[0]

#convert years to months
kydat[:,:,0]=(kydat[:,:,0]-2005)*12
jydat[:,:,0]=(jydat[:,:,0]-2005)*12

#%%
#function required for wst test
def delta(data,error):
    mean = np.sum(data*error**-2)/np.sum(error**-2)
    return (data-mean)/error

#wst function
def wstest(jdata,jerr,kdata,kerr):
    kdelta = delta(kdata,kerr)
    jdelta = delta(jdata,jerr)
    return (1/np.sqrt(7*6))*np.sum(kdelta*jdelta)

#Returns the number of points adjacent to it within 3stds will be 1 if none match becase its matching itself
def wstestremove(jdata,jerr,kdata,kerr):
    wst = np.zeros((114243,7))
    for j in range(7):
        currj = np.delete(jdata,j)
        currk = np.delete(kdata,j)
        currjerr = np.delete(jerr,j)
        currkerr = np.delete(kerr,j)
        print(j)
        for i in range(114243):
            wst[i,j]=wstest(currj,currjerr,currk,currkerr)
    return wst

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

def lightcurve(id):
    fig,ax = plt.subplots(2)
    fig.tight_layout()
    ax[1].plot(kydat[id,:,0],kydat[id,:,1],'b-')
    ax[1].errorbar(kdat[id,:,0],kdat[id,:,1],yerr=kdat[id,:,2],color="red",marker="x",lw=0,elinewidth=1,capsize=1.5)
    ax[0].plot(jydat[id,:,0],jydat[id,:,1],'b-')
    ax[0].errorbar(jdat[id,:,0],jdat[id,:,1],yerr=jdat[id,:,2],color="red",marker="x",lw=0,elinewidth=1,capsize=1.5)
    ax[1].set(ylabel="K")
    ax[0].set(ylabel="J")
    ax[0].set(xlabel="months")
    ax[0].plot(jydat[id,:,0],np.ones(7)*np.median(jdat[id,:,1]),'g-')
    ax[1].plot(kydat[id,:,0],np.ones(7)*np.median(kdat[id,:,1]),'g-')
    fig.suptitle(str(("DR11= " , int(lut[id,0]))))

#%%
#peak value dector combined with wst to return variable objects with single peaks

#parameters:
wstLimit = 0
adjLimit = 1

"""#monthy wst
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

objType = np.zeros(114243)

counts=[0,0,0]#sn/agn count

for i in range(114243):
    #if object is variable
    if wstest(jydat[i,:,1],jydat[i,:,2],kydat[i,:,1],kydat[i,:,2])>wstLimit:
        jadj,kadj = adjacencytest(jydat[i,:,1],jydat[i,:,2],kydat[i,:,1],kydat[i,:,2],adjLimit)
        counts[0]+=1
        if np.any(arrAGN==i):
            counts[2]+=1
            objType[i] = 2
        elif np.sum(jadj)==1 and np.sum(kadj)==1:
            objType = 1
            counts[1]+=1
            lightcurve(i)
            plt.show()
            print(i)


print("total,SN, AGN = ", counts)

"""
import matplotlib.pyplot as plt
for id in range(np.size(kpeakObjs[:,0,0])):
   
    plt.show()
"""

#%%
#functions for stat test analysis.

#returns the number of stds a point is from the mean
def peaktest(jdata,kdata,jerr,kerr):
    javg = np.mean(jdata)
    kavg = np.mean(kdata)
    
    jpeaks = (jdata-javg)/jerr
    kpeaks = (kdata-kavg)/kerr

    jdeviations = np.abs(jpeaks)
    kdeviations = np.abs(kpeaks)

    return jpeaks,kpeaks,jdeviations,kdeviations


#wst but removing a point YEARLY DATA BEST



#%%
i = 95788
lightcurve(i)
print(kdiffnum[i],jdiffnum[i])
print(objType[i])
print(lut[i,0])

plt.show()


#%%

chisq = wst
import matplotlib.pyplot as plt
plt.hist((chisq))
plt.show()
plt.plot(chisq)
plt.show()

std = np.std(chisq)
mean = np.mean(chisq)

big = np.where(np.logical_and(chisq>3*std+mean,True))#,chisq<10000000*std+mean))

kdatbig = kdat[big]
jdatbig =jdat[big]

chisqlst = chisq[big]
indx = np.arange(kdat[:,0,0].size)[big]

for i in range(kdatbig[:,0,0].size):
    plt.plot(kdatbig[i,:,0],kdatbig[i,:,1],'b.')
    plt.errorbar(kdatbig[i,:,0],kdatbig[i,:,1],kdatbig[i,:,2],color='blue',linestyle='')
    
    plt.plot(jdatbig[i,:,0],jdatbig[i,:,1],'r.')
    plt.errorbar(jdatbig[i,:,0],jdatbig[i,:,1],jdatbig[i,:,2],color='red',linestyle='')

    print((chisqlst[i]-mean)/std,indx[i])
    plt.show()
lut = np.loadtxt("data/LUT.npy")
ids = lut[big]



#%%
#plotting of spatial region of "zero" objects
kzero = np.any(kdat[:,:,1]==0,axis=1)
jzero = np.any(jdat[:,:,1]==0,axis=1)

posk=lut[kzero,:]
posj=lut[jzero,:]


def setMask(array): #Funcion to remove values outside of the range specified
    array= array.T
    xmin,ymin = 34.2,-5.4
    xmax,ymax = 34.8,-4.8
    mask = np.logical_or(array[:,0]<xmin,array[:,0]>xmax)
    mask |= np.logical_or(array[:,1]<ymin,array[:,1]>ymax)
    retarray = array[~mask]
    return retarray.T

xmin,ymin = 34.2,-5.4
xmax,ymax = 34.8,-4.8

kcoremask = np.logical_or(posk[:,1]<xmin,posk[:,1]>xmax)
kcoremask |= np.logical_or(posk[:,2]<ymin,posk[:,2]>ymax)

jcoremask = np.logical_or(posj[:,1]<xmin,posj[:,1]>xmax)
jcoremask |= np.logical_or(posj[:,2]<ymin,posj[:,2]>ymax)



corek=posk[~kcoremask]
corej=posj[~jcoremask]

plt.figure(dpi=120)
plt.plot(lut[:,1],lut[:,2],'b,')
plt.plot(posk[:,1],posk[:,2],'r.')
plt.plot(posj[:,1],posj[:,2],'r.')
plt.plot(corek[:,1],corek[:,2],'g.')
plt.plot(corej[:,1],corej[:,2],'g.')
plt.show()
#%%
#lightcurve of the field and variable objects within it.
j=0
tot = np.zeros(7)
toterr = tot.copy()
for i in range(114243):
    if wstest(jydat[i,:,1],jydat[i,:,2],kydat[i,:,1],kydat[i,:,2])>1:
        toterr += kydat[i,:,2]
        tot += kydat[i,:,1]
        j=j+1
print(j)
plt.errorbar(jydat[0,:,0],tot,toterr)


#%%
print(np.where(lut==62253))
lightcurve(95788)


#%%
#Varability imbalance detection
"""
kmean = np.expand_dims(np.mean(kdat[:,:,1],axis=1),axis=-1)
jmean = np.expand_dims(np.mean(jdat[:,:,1],axis=1),axis=-1)

kUpOrDown = (kdat[:,:,1] > kmean) #create a masking array to apply to the original data to seperate into >mean and <mean
jUpOrDown = (jdat[:,:,1] > jmean)

print(np.sum(kUpOrDown)/np.size(kdat))
#this is a really stupid way of doing it but its funciton and doesnt work so I cba to change it
kdist = np.zeros(114243)
jdist = np.zeros(114243)
for i in range(114243):
    kup   = kdat[i,kUpOrDown[i,:],1]
    jup   = jdat[i,jUpOrDown[i,:],1]
    kdown = kdat[i,~kUpOrDown[i,:],1]
    jdown = jdat[i,~jUpOrDown[i,:],1]

    kuperr   = kdat[i,kUpOrDown[i,:],2]
    juperr   = jdat[i,jUpOrDown[i,:],2]
    kdownerr = kdat[i,~kUpOrDown[i,:],2]
    jdownerr = jdat[i,~jUpOrDown[i,:],2]

    kdist[i] = kdist[i] + (np.sum((kmean[i]-kup)/kuperr) + np.sum((kmean[i]-kdown)/kdownerr))/1#kmean[i]
    jdist[i] = jdist[i] + (np.sum((jmean[i]-jup)/juperr) + np.sum((jmean[i]-jdown)/jdownerr))/1#jmean[i]

plt.hist(np.log(kdist),bins=1000,alpha=0.6)
plt.hist(np.log(jdist),bins=1000,alpha=0.6)
plt.show()

for i in range(114243):
    if kdist[i]>10:
        plt.plot(kydat[i,:,0],kydat[i,:,1],'r-')
        plt.errorbar(kdat[i,:,0],kdat[i,:,1],yerr=kdat[i,:,2],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
        print(i,kdist[i])
        plt.show()
"""



#%%

for i in arrSN:
    lightcurve(i)
    plt.show()

#%%
    
lightcurve(95788)