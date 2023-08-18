import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util_os import util_os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

import glob

# flux density
def getFluxDensity(hi, dv, D):
    S21 = []
    for m in hi:
        s = (1/(2.36e5)) * (m/(1)) * (1/dv)*((1/D)**2)
        S21.append(s)
    return np.array(S21)

def makeHIprofile(gal):

    cdmPath = util_os.getfilepath_cdm(gal, 0)
    hifiles = util_os.getFilePath_HI(gal)

    # flux density params
    D = 70 #Mpc (distance observed at; note this is kinda arbitrary atm)
    dv= 11.2 #km/s (velocity resolution)

    # style params
    tsize = 24
    asize = 18

    lw = 5  # profile linewidth
    lwW = 2 # width linewidth

    cprof = 'royalblue'
    cwid = 'magenta'

    # widths to measure
    Ws = [50,20,10]

    fcdm = fits.open(hifiles[1], do_not_scale_image_data=True)

    bscale = fcdm[0].header['BSCALE']
    bzero = fcdm[0].header['BZERO']
    blank = fcdm[0].header['BLANK']
    data = fcdm[0].data

    vstart = fcdm[0].header['CRVAL3']
    dv = fcdm[0].header['CDELT3']

    # initial scaling
    data = data * bscale + bzero
    data = np.power(np.full(data.shape,10), data)
    data[np.where(abs(data - 10**(blank*bscale + bzero))<10e-5)]=0.0

    # convert from solar mass to Jy (el bhadri et al 2018 for mock obs at 70 MPC)
    S = (1/(2.36e5))*(data)*(1/dv)*(D**-2)

    # velocities
    xaxis = np.arange(vstart, vstart+(dv*len(S[:,0,0])), dv)

    # sum along spectral axis
    hiprof = np.sum(np.sum(S,1),1)

    plt.plot(xaxis, hiprof, c=cprof, linewidth=4)

    # do line widths of
    vmax = max(hiprof)
    wids = [-1,-1,-1,-1,-1,-1]
    Ws = [10,20,50]
    for j,p in enumerate(Ws):
        val = (p/100)*vmax

        idxs_mass = np.argwhere(np.diff(np.sign(hiprof - np.full(len(hiprof),
            val)))).flatten()
        if len(idxs_mass) < 2: continue
        x1 = (xaxis[idxs_mass[0]]+xaxis[idxs_mass[0]+1])/2
        x2 = (xaxis[idxs_mass[1]]+xaxis[idxs_mass[1]+1])/2

        width = x2-x1
        wids[j] = width
        #print('CDM W'+str(p)+'_'+orientations[i]+': '+str(width))

        plt.plot([x1,x2],[val,val], linewidth=3, zorder=2, c=cwid,
            label='_nolegend_')
        plt.scatter([x1,x2],[val,val], s=50, label='_nolegend_', zorder=2,
            marker='o', c=cwid)

    # prettify
    plt.title('r'+str(gal), fontfamily='serif', fontsize=24, weight='bold', loc='left')
    plt.xlabel('Gas Velocity [km/s]', fontsize=20)
    plt.ylabel('Flux Density [Jy]', fontsize=20)

    plt.savefig('../figures/HIProfiles/appendix/r'+str(gal)+'.png', dpi=600)

# end makeHIprofile

# get haloIDs
cdmHalos = util_os.getGalaxies()[0]

print('Making HI Profiles')
# make figs
f=open('/home/jw1624/H1-merian/csvs/HI_widths.txt', 'w')
f.write('galaxy,w50_cdm,w20_cdm,w10_cdm,')
f.write('w50_sidm,w20_sidm,w10_sidm,')
f.write('dW_cdm,dW_sidm,')
f.write('K_cdm,K_sidm\n')
f.close()
for g in cdmHalos:
    print(' halo '+str(g)+'...', end='')
    makeHIprofile(g)
    print('done')
