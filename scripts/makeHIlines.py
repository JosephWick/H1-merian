import sys
sys.path.insert(0, '../util/')
from util import util

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import glob

def makeHIprofile(hID, withSIDM=False):
    dpath = '/home/jw1624/H1-merian/h1lines/'

    # style params
    tsize = 24
    asize = 18

    lw = 5
    lwW = 2

    cCDM = 'tab:blue'
    cCDMw = 'b'
    cSIDM = 'tab:orange'
    cSIDMw = 'r'

    orientations = ['x','y','z']

    # widths to measure
    Ws = [50,20,10]

    # fig setup
    fig,axs = plt.subplots(1,3, figsize=(14,5))
    plt.figure(facecolor='white')

    # get the data
    fcdm = glob.glob(dpath+'r'+str(hID)+'_cdm*')
    fsidm= glob.glob(dpath+'r'+str(hID)+'_sidm*')

    fcdm.sort()
    fsidm.sort()

    # plot
    plt.suptitle('HI Profile for Galaxy '+str(hID), fontsize=tsize)
    for i in range(len(fcdm)):
        cdmx = pd.read_csv(fcdm[i], sep='\s+', header=None)[0]
        cdmy = pd.read_csv(fcdm[i], sep='\s+', header=None)[1]

        sidmx = pd.read_csv(fsidm[i], sep='\s+', header=None)[0]
        sidmy = pd.read_csv(fsidm[i], sep='\s+', header=None)[1]

        axs[i].plot(cdmx, cdmy, linewidth=lw+1, c=cCDM)
        axs[i].plot(sidmx, sidmy, linewidth=lw-1, c=cSIDM)

        axs[i].set_title(orientations[i]+' Orientation')
        axs[i].set_xlabel(r'Velocity [km/s]', fontsize=asize)
        axs[i].set_ylabel(r'Mass [$M_\odot$]', fontsize=asize)

        # do line widths of cmd
        vmax = max(cdmy)
        for j,p in enumerate(Ws):
            val = (p/100)*vmax

            idxs_mass = np.argwhere(np.diff(np.sign(cdmy - np.full(len(cdmy),val)))).flatten()
            if len(idxs_mass) < 2: continue
            x1 = (cdmx[idxs_mass[0]]+cdmx[idxs_mass[0]+1])/2
            x2 = (cdmx[idxs_mass[1]]+cdmx[idxs_mass[1]+1])/2

            width = x2-x1
            print('CDM W'+str(p)+'_'+orientations[i]+': '+str(width))

            axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cCDMw)
            axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2, c=cCDMw)

        # do line widths of sidm
        if withSIDM:
            vmax = max(sidmy)
            for j,p in enumerate(Ws):
                val = (p/100)*vmax

                idxs_mass = np.argwhere(np.diff(np.sign(sidmy - np.full(len(sidmy),val)))).flatten()
                if len(idxs_mass) < 2: continue
                x1 = (sidmx[idxs_mass[0]]+sidmx[idxs_mass[0]+1])/2
                x2 = (sidmx[idxs_mass[1]]+sidmx[idxs_mass[1]+1])/2

                width = x2-x1
                print('SIDM W'+str(p)+'_'+orientations[i]+': '+str(width))

                axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cSIDMw)
                axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2, marker=markers[j], c=cSIDMw)

    plt.tight_layout()
    if withSIDM:
        plt.savefig('../figures/HIProfiles/HI_'+str(hID)+'_2.png')
    else:
        plt.savefig('../figures/HIProfiles/HI_'+str(hID)+'.png')
# end makeHIprofile

# get haloIDs
cdmHalos, sidmHalos = util.getGalaxies()

print('Making HI Profiles')
# make figs
for g in cdmHalos:
    print(' halo '+str(g)+'...', end='')
    if g in sidmHalos:
        makeHIprofile(g, withSIDM=True)
    makeHIprofile(g)
    print('done')
