#%%
import numpy as np
import matplotlib.pyplot as plt
#import dataset

kdat = np.loadtxt("data/kdata.npy").reshape((114243,38,3))
jdat = np.loadtxt("data/jdata.npy").reshape((114243,35,3))

kydat = np.loadtxt("data/kydata.npy").reshape((114243,7,3))
jydat = np.delete(np.loadtxt("data/jydata.npy").reshape((114243,8,3)),(1),axis=1)

#import lookup tables (arrays of DR11 ids)
lut = np.loadtxt("data/LUT.npy")
AGNlut = np.loadtxt("data/AGNLUT.npy")

#find suspected SN array indexes and DR11 ids
foundSN = np.array([147703,89852,62253,66612,24605,219616,148305,263613,214136,210486,154039,240562,224386,217727,117014,252582,248061])
arrSN=foundSN
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

#Define tests
def wstest(jdata,kdata,jerr,kerr):
    size = len(jdata)
    kdelta = (kdata-np.sum(kdata*kerr**-2)/np.sum(kerr**-2))/kerr
    jdelta = (jdata-np.sum(jdata*jerr**-2)/np.sum(jerr**-2))/jerr
    index = (1/np.sqrt(size*(size-1)))*np.sum(kdelta*jdelta)
    return index


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
        r=True
    else: r=False

    #if total error is greater below mean than above then test is positive. 
    if np.sum(jpeaks)<-1 and np.sum(kpeaks)<-1:
        s=True
    else:s=False

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
def wstestremove(jdata,kdata,jerr,kerr):
    size = np.size(jdata)
    wst = np.zeros(size)
    for i in range(size):
        currj = np.delete(jdata,i)
        currk = np.delete(kdata,i)
        currjerr = np.delete(jerr,i)
        currkerr = np.delete(kerr,i)
        wst[i]=wstest(currj,currk,currjerr,currkerr)
    #plt.plot(wst)
    if np.sum(wst<0.7)==1 and np.std(wst)>1: 
        r = 1
    else: r=0
    return r, np.std(wst)

jmag = 30-2.5*np.log(jdat)
kmag = 30-2.5*np.log(kdat)
jymag = 30-2.5*np.log(jydat)
kymag = 30-2.5*np.log(kydat)

#magnitude lightcurve plotter
def lightcurve2(id):
    
    fig,ax = plt.subplots(2)
    plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.94,hspace=0)
    
    ax[1].plot(jydat[id,:,0],kymag[id,:,1],'b-')
    ax[0].plot(jydat[id,:,0],jymag[id,:,1],'b-')
    
    #error bars
    ax[0].errorbar(jdat[id,:,0] ,jmag[id,:,1] ,yerr=(jdat[id,:,2] /jdat [id,:,1])*jmag [id,:,1],color="red" ,marker="x",lw=0,markersize =4,elinewidth=0.5,capsize=0)
    ax[0].errorbar(jydat[id,:,0],jymag[id,:,1],yerr=(jydat[id,:,2]/jydat[id,:,1])*jymag[id,:,1],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
    ax[1].errorbar(kdat[id,:,0] ,kmag[id,:,1] ,yerr=(kdat[id,:,2] /kdat [id,:,1])*kmag [id,:,1],color="red" ,marker="x",lw=0,markersize=4,elinewidth=0.5,capsize=0)
    ax[1].errorbar(kydat[id,:,0],kymag[id,:,1],yerr=(kydat[id,:,2]/kydat[id,:,1])*kymag[id,:,1],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
        
    #titles
    ax[0].set(ylabel="J")
    ax[1].set(ylabel="K")
    ax[1].set(xlabel="Month")
    ax[0].xaxis.set_visible(False)
    fig.suptitle(("DR11= %i" %(int(lut[id,0]))))
    
    #scale timeline
    ax[0].set_xlim(-2,87)
    ax[1].set_xlim(-2,87)

    #mean month line
    ax[0].plot(jydat[id,:,0],np.ones(7)*np.mean(jmag[id,:,1]),'g-')
    ax[1].plot(kydat[id,:,0],np.ones(7)*np.mean(kmag[id,:,1]),'g-')

    ax[0].invert_yaxis()
    ax[1].invert_yaxis()


def lightcurve(id):
    
    fig,ax = plt.subplots(2)
    plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.94,hspace=0)
    
    ax[1].plot(jydat[id,:,0],kydat[id,:,1],'b-')
    ax[0].plot(jydat[id,:,0],jydat[id,:,1],'b-')
    
    #error bars
    ax[0].errorbar(jdat[id,:,0] ,jdat[id,:,1] ,yerr=jdat[id,:,2] ,color="red",marker="x",markersize =4,lw=0,elinewidth=0.5,capsize=0)
    ax[0].errorbar(jydat[id,:,0],jydat[id,:,1],yerr=jydat[id,:,2],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
    ax[1].errorbar(kdat[id,:,0] ,kdat[id,:,1] ,yerr=kdat[id,:,2] ,color="red",marker="x",markersize=4,lw=0,elinewidth=0.5,capsize=0)
    ax[1].errorbar(kydat[id,:,0],kydat[id,:,1],yerr=kydat[id,:,2],color="blue",marker="x",lw=0,elinewidth=1,capsize=1.5)
        
    #titles
    ax[0].set(ylabel="J")
    ax[1].set(ylabel="K")
    ax[1].set(xlabel="Month")
    ax[0].xaxis.set_visible(False)
    fig.suptitle(("DR11= %i" %(int(lut[id,0]))))
    
    #scale timeline
    ax[0].set_xlim(-2,87)
    ax[1].set_xlim(-2,87)

    #mean month line
    ax[0].plot(jydat[id,:,0],np.ones(7)*np.mean(jydat[id,:,1]),'g-')
    ax[1].plot(kydat[id,:,0],np.ones(7)*np.mean(kydat[id,:,1]),'g-')

#%%
#peak value dector combined with wst to return variable objects with single peaks

#0 = peaktest, 1 = balance test, 
counts=np.zeros((5))#sn/agn count
wststd = np.array([])
for i in range(114243):
    #if object is variable
    countprev = counts.copy()
    if wstest(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])>-1000:
        #test each point to:
        r,s =peaktest(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])
        counts[0] += r
        counts[1] += s
        counts[2] += adjacencytest(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])
        t, std = wstestremove(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])
        counts[3]+=t
        if t==1:
            wststd=np.append(wststd,std)
        if  np.sum(counts[:-1]-countprev[:-1])>0:
            counts[4] += 1
        if counts[3]-countprev[3] ==1:
            lightcurve(i)
            plt.show()
            print(int(lut[i,0]))
            print(counts[:-1]-countprev[:-1])

print(counts)

#%%
#SN shortlist
peakTestSN = np.array([148305, 145425, 117014, 97639, 62253, 13145, 66612, 24605, 219616, 214136, 292377, 202060, 153802, 147703])
balanceTestSN = np.array([148305, 145425, 62253])
adjTestSN = np.array([148305, 117014, 62253, 13145, 147703, 66612, 24065, 154039, 202060, 153802])
wstRemoveSn = np.array([104987,62253,13145,147703,66612,292377,217727,202060,153802])#([62253, 13145, 147703, 66612, 219616, 252582, 202060, 153802])

#153802 111791 62253 104987 147703

#                                                             
#                                                             weird                              
#            Certainty?:   0      0       1      0     0.5     0     0       0        1      1        1       0.2     ?       0       1       1       1   
SNShortlist = np.array([148305, 117014, 62253, 13145, 66612, 24605, 219616, 202060, 153802, 147703, 214136, 292377, 252582, 154039, 104987, 111791, 210486 ])
SNXrayDet   = np.array([True  , True  , False, False, False, False, False , True  , False , False , False , False , False , False , False , False , False ])
SNpz        = np.array([               1.5156,0.7362,1.0109,1.9400,0.9457         , 2.7708, 1.5895, 1.4412, 1.2946, 0.9206, 3.0199, 0.2993, 0.1444, 2.3681])
arrSN = SNShortlist.copy()

SNconfident = SNShortlist[SNXrayDet]

for i in range(np.size(SNShortlist)):
    if SNXrayDet[i] == True:
        print(SNShortlist[i])
        print(SNXrayDet[i])

        a, = np.where(lut[:,0].astype(int)==SNShortlist[i])[0]
        arrSN[i] = a
        lightcurve(a)
        plt.show()
print(np.sum(~SNXrayDet))

"""
import matplotlib.pyplot as plt
for id in range(np.size(kpeakObjs[:,0,0])):
   
    plt.show()
"""

#%% known object
#i = 95788 #found sn
#i = 40030
i = np.where(lut==194551)[0][0]
wstestremove(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])
#r,wst = wstestremove(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])

lightcurve(i)
#lightcurve2(i)
#plt.show()

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
plt.title("Objects with 0 values")
plt.xlabel("RA/deg")
plt.ylabel("Dec/deg")
plt.show()
#%%
#lightcurve of the field and variable objects within it.
j=0
tot = np.zeros(7)
toterr = tot.copy()
for i in range(114243):
    if wstest(jydat[i,:,1],kydat[i,:,1],jydat[i,:,2],kydat[i,:,2])>0.7:
        toterr += kydat[i,:,2]
        tot += kydat[i,:,1]
        j=j+1
print(j)
plt.errorbar(jydat[0,:,0],tot,toterr)


#%%

superColour = np.loadtxt("data/superColour.npy")

#print(np.where(superColour[:,0]==lut[:,0])[0])



#convert agn dr11s to array indexes
SCarr = SNconfident.copy()
for i in range(np.size(SCarr)):
    SCarr[i]=np.where(SNconfident[i]==superColour[:,0])[0].astype(int)

for i in SCarr:
    print(superColour[i,:].astype(int))

#8xSF (2) 1xpassive (1), 2xundef (neg)
