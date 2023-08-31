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

# mu_n for kurtosis (defn from El Bhadri et al 2018)
def mu_n(v,s, n):
    # s is line profile
    # v is velocity bins
    vbar = sum(np.array(s*v))/sum(s)

    top = (s*((v-vbar)**n)).sum()
    bottom = s.sum()
    return top/bottom

# getKurtosis
def getKurtosis(v,s):
    return mu_n(v,s,4)/(mu_n(v,s,2)**2) - 3

def makeHIprofile(gal, withSIDM=False, doExport=True):
    f=open('/home/jw1624/H1-merian/csvs/HI_widths.txt', 'a')

    h1files_cdm = util_os.getFilePath_HI_cdm(gal)
    h1files_sidm = util_os.getFilePath_HI_sidm(gal)

    # flux density params
    D = 70 #Mpc (distance observed at; note this is kinda arbitrary atm)
    dv= 11.2 #km/s (velocity resolution)

    # style params
    tsize = 24
    asize = 18

    lw = 5  # profile linewidth
    lwW = 2 # width linewidth

    cCDM = 'firebrick'
    cCDMw = 'r'
    cSIDM = 'royalblue'
    cSIDMw = 'b'

    orientations = ['x','y','z']
    markers=['o','^','x']

    # widths to measure
    Ws = [50,20,10]

    # fig setup
    fig,axs = plt.subplots(1,3, figsize=(14,5), facecolor='white')

    # plot
    plt.suptitle('HI Profile for Galaxy '+str(gal), fontsize=tsize)
    K_cdm = -1
    for i in range(len(h1files_cdm)):
        # read data
        # no scaling so I can follow Alyson's IDL code
        fcdm = fits.open(h1files_cdm[i], do_not_scale_image_data=True)
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

        K_cdm = getKurtosis(xaxis, hiprof)

        axs[i].plot(xaxis, hiprof, linewidth=lw+1, c=cCDM)

        axs[i].set_title(orientations[i]+' Orientation')
        axs[i].set_xlabel(r'Velocity [km/s]', fontsize=asize)
        axs[i].set_ylabel('Flux Density [Jy]', fontsize=asize)

        # do line widths of cmd
        vmax = max(hiprof)
        wids = [-1,-1,-1,-1,-1,-1]
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

            axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cCDMw,
                label='_nolegend_')
            axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2,
                marker=markers[j], c=cCDMw)

        # do line widths of sidm
        K_sidm = -1
        if withSIDM:
            # read data
            # no scaling so I can follow Alyson's IDL code precisely
            fsidm = fits.open(h1files_sidm[i], do_not_scale_image_data=True)
            bscale = fsidm[0].header['BSCALE']
            bzero = fsidm[0].header['BZERO']
            blank = fsidm[0].header['BLANK']
            data = fsidm[0].data

            vstart = fsidm[0].header['CRVAL3']
            dv = fsidm[0].header['CDELT3']

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

            K_sidm = getKurtosis(xaxis, hiprof)

            axs[i].plot(xaxis, hiprof, linewidth=lw-1, c=cSIDM)

            vmax = max(hiprof)
            for j,p in enumerate(Ws):
                val = (p/100)*vmax

                idxs_mass = np.argwhere(np.diff(np.sign(hiprof - np.full(len(hiprof),
                    val)))).flatten()
                if len(idxs_mass) < 2: continue
                x1 = (xaxis[idxs_mass[0]]+xaxis[idxs_mass[0]+1])/2
                x2 = (xaxis[idxs_mass[1]]+xaxis[idxs_mass[1]+1])/2

                width = x2-x1
                wids[j+3] = width
                #print('SIDM W'+str(p)+'_'+orientations[i]+': '+str(width))

                axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cSIDMw,
                    label='_nolegend_')
                axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_',
                    zorder=2, marker=markers[j], c=cSIDMw)

        if doExport:
            f.write(str(gal)+'_'+orientations[i]+',')
            for w in range(len(wids)-1):
                f.write(str(wids[w])+',')
            f.write(str(wids[-1])+',')

            # calculate dW
            # dW = (W_20 - W_50)/W_50 (el Bhadri 2018)
            dW_cdm = -1
            dW_sidm= -1
            if wids[0] != -1 and wids[1] != -1:
                dW_cdm = (wids[1] - wids[0])/wids[0]
            if wids[3] != -1 and wids[4] != -1:
                dW_sidm = (wids[4] - wids[3])/wids[3]
            f.write(str(dW_cdm)+','+str(dW_sidm)+',')


    f.close()
    plt.tight_layout()
    if withSIDM:
        axs[0].legend(['CDM','SIDM'])
        plt.savefig('/home/jw1624/H1-merian/figures/HIProfiles/HI_'+str(gal)+'_2.png')
    else:
        plt.savefig('/home/jw1624/H1-merian/figures/HIProfiles/HI_'+str(gal)+'.png')

# end makeHIprofile

# get haloIDs
cdmHalos, sidmHalos, adiabaticHalos = util_os.getGalaxies()

print('Making HI Profiles')
# make figs
f=open('/home/jw1624/H1-merian/csvs/HI_widths.txt', 'w')
f.write('galaxy,w50_cdm,w20_cdm,w10_cdm,')
f.write('w50_sidm,w20_sidm,w10_sidm,')
f.write('dW_cdm,dW_sidm,')
f.close()
for g in cdmHalos:
    print(' halo '+str(g)+'...', end='')
    if g in sidmHalos:
        makeHIprofile(g, withSIDM=True, doExport=True)
        makeHIprofile(g, withSIDM=False,doExport=False)
    else:
        makeHIprofile(g, withSIDM=False,doExport=True)
    print('done')
