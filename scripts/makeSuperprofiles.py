import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import numpy as np
import pandas as pd

from astropy.io import fits

import matplotlib.pyplot as plt

import glob
import math

import scipy.optimize

# for fitting HI profile
def gaussHermite(x, A, s, h, mu):
    y = (x-mu)/s
    outside = A*np.exp((-1*np.power(y,2))/2)
    inside = 1 + (h/np.sqrt(6))*(2*np.sqrt(2)*np.power(y,3) - 3*np.sqrt(2)*y)
    return outside*inside

# for fitting superprofile
def gaussian(x, A, mu, s):
    return (A/(s*np.sqrt(2*np.pi))) * np.exp((-1*np.power(x-mu,2))/(2*np.power(s,2)))

# makes figures
def superprofile(hID, withSIDM=False):
    f=open('/home/jw1624/H1-merian/csvs/superProfiles.txt', 'a')

    cdmPath, sidmPath, _ = util.getfilepath(hID)

    h1files = glob.glob(cdmPath+'/*.fits')
    h1files.sort()
    if len(h1files) == 4: h1files = h1files[1:]

    # load data (yhat orientation for now)
    f = fits.open(h1files[1], do_not_scale_image_data=True)
    bscale = f[0].header['BSCALE']
    bzero = f[0].header['BZERO']
    blank = f[0].header['BLANK']
    data = f[0].data

    vstart = f[0].header['CRVAL3']
    dv = f[0].header['CDELT3']
    D = 70 #mpc

    # initial scaling
    data = data * bscale + bzero
    data = np.power(np.full(data.shape,10), data)
    data[np.where(abs(data - 10**(blank*bscale + bzero))<10e-5)]=0.0

    # convert from solar mass to Jy (el bhadri et al 2018 for mock obs at 70 MPC)
    S = (1/(2.36e5))*(data)*(1/dv)*(D**-2)

    for i in range(len(S[0,:,0])):
        for j in range(len(S[0,0,:])):
            xaxis = np.arange(vstart, vstart+(dv*len(S[:,i,j])), dv)
            # get naive vmax
            vmax = xaxis[np.where(S[:,i,j] == S[:,i,j].max())][0]

            # do fit
            if max(S[:,i,j]) > 10e-6: #10
                maxval = S[:,i,j].max()
                p0 = [maxval, 10, 0, xaxis[np.where(S[:,i,j]==maxval)][0]]
                p,_ = scipy.optimize.curve_fit(gaussHermite, xaxis, S[:,i,j], maxfev=100000, p0=p0)
                fity = gaussHermite(xaxis, p[0],p[1],p[2],p[3])
                vmaxfit = xaxis[np.where(fity == fity.max())][0]

                # shift array such that vmax is over zero
                S[:,i,j] = np.roll(S[:,i,j], math.ceil(-vmaxfit/dv))

    # sum along spectral axis
    hiprof = np.sum(np.sum(S,1),1)

    # fit gaussian
    p0 = [max(hiprof), 0, 100]
    p,_ = scipy.optimize.curve_fit(gaussian, xaxis, hiprof, maxfev=100000, p0=p0)
    g = gaussian(xaxis, p[0],p[1],p[2])

    # get superprofile params
    sigmaCentral = p[1]
    hwhm = 2*np.sqrt(2*np.log(2))*sigmaCentral

    fwings = 0
    sumS = 0
    for i,v in enumerate(xaxis):
        if abs(v)>hwhm:
            fwings += (S[i] - g(v, p[0],p[1],p[2]))
            sumS += S[i]
    fwings = fwings/sumS

    sigmaWings = 0
    sumSmG = 0
    for i,v in enumerate(xaxis):
        if abs(v) > hwhm:
            sigmaWings += (S[i] - g(v,p[0],p[1],p[2]))*(v**2)
            sumSmG += S[i] - g(v,p[0],p[1],p[2])
    sigmaWings = np.sqrt(sigmaWings/sumSmG)

    a = 0
    for i,v in enumerate(xaxis):
        if abs(v) > hwhm:
            a += np.sqrt((S[i] - S[-i-1])**2)
    a = a/sumSmG

    # write to table
    f.write(str(hID)+','+str(sigmaCentral)+','+str(fwings)+','
        +str(sigmaWings)+','+str(a)+'\n')
    f.close()

    # do figure
    fig = plt.figure(figsize=(8,8), facecolor='w')

    plt.plot(xaxis, hiprof, linewidth=2, c='k')
    plt.plot(xaxis, g, linewidth=2, linestyle='r--')
    plt.fill_between(xaxis, hiprof, g, color='g', alpha=0.5)

    # prettify
    plt.title('Halo '+str(hID), fontsize=22)
    plt.xlabel('Offset Velocity [km/s]', fontsize=20)
    plt.ylabel('Flux [Jy]', fontsize=20)
    plt.legend(['Gaussian Fit', 'HI Superprofile'])

    plt.savefig('/home/jw1624/H1-merian/figures/superprofiles/r'+str(hID)+'.png')

    # end



# get galaxies
cdmHalos,sidmHalos,_ = util.getGalaxies()

# set up file
f=open('/home/jw1624/H1-merian/csvs/superProfiles.txt', 'w')
f.write('galaxy,sigmaCentral,fWings,sigmaWings,a')
f.close()

for g in cdmHalos:
    if g in sidmHalos:
        superprofile(g, withSIDM=True)
    superprofile(g, withSIDM=False)
    print('halo '+str(g)+'... done')
