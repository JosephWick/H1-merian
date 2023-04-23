import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import glob

def makeHIprofile(hID, withSIDM=False):
    f=open('/home/jw1624/H1-merian/h1lines/widths.txt', 'a')

    dpath = '/home/jw1624/H1-merian/h1lines/'

    # style params
    tsize = 24
    asize = 18

    lw = 5
    lwW = 2

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

    # get the data
    fcdm = glob.glob(dpath+'r'+str(hID)+'_cdm*')
    fsidm= glob.glob(dpath+'r'+str(hID)+'_sidm*')

    fcdm.sort()
    fsidm.sort()

    # plot
    plt.suptitle('HI Profile for Galaxy '+str(hID), fontsize=tsize)
    for i in range(len(fcdm)):
        f.write(str(hID)+'_'+orientations[i]+'\t')

        cdmx = pd.read_csv(fcdm[i], sep='\s+', header=None)[0]
        cdmy = pd.read_csv(fcdm[i], sep='\s+', header=None)[1]

        axs[i].plot(cdmx, cdmy, linewidth=lw+1, c=cCDM)

        axs[i].set_title(orientations[i]+' Orientation')
        axs[i].set_xlabel(r'Velocity [km/s]', fontsize=asize)
        axs[i].set_ylabel(r'Mass [$M_\odot$]', fontsize=asize)

        # do line widths of cmd
        vmax = max(cdmy)
        wids = [-1,-1,-1,-1,-1,-1]
        for j,p in enumerate(Ws):
            val = (p/100)*vmax

            idxs_mass = np.argwhere(np.diff(np.sign(cdmy - np.full(len(cdmy),val)))).flatten()
            if len(idxs_mass) < 2: continue
            x1 = (cdmx[idxs_mass[0]]+cdmx[idxs_mass[0]+1])/2
            x2 = (cdmx[idxs_mass[1]]+cdmx[idxs_mass[1]+1])/2

            width = x2-x1
            wids[j] = width
            #print('CDM W'+str(p)+'_'+orientations[i]+': '+str(width))

            axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cCDMw)
            axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2, marker=markers[j], c=cCDMw)

        # do line widths of sidm
        if withSIDM:
            sidmx = pd.read_csv(fsidm[i], sep='\s+', header=None)[0]
            sidmy = pd.read_csv(fsidm[i], sep='\s+', header=None)[1]

            axs[i].plot(sidmx, sidmy, linewidth=lw-1, c=cSIDM)

            vmax = max(sidmy)
            for j,p in enumerate(Ws):
                val = (p/100)*vmax

                idxs_mass = np.argwhere(np.diff(np.sign(sidmy - np.full(len(sidmy),val)))).flatten()
                if len(idxs_mass) < 2: continue
                x1 = (sidmx[idxs_mass[0]]+sidmx[idxs_mass[0]+1])/2
                x2 = (sidmx[idxs_mass[1]]+sidmx[idxs_mass[1]+1])/2

                width = x2-x1
                wids[j+3] = width
                #print('SIDM W'+str(p)+'_'+orientations[i]+': '+str(width))

                axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cSIDMw)
                axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2, marker=markers[j], c=cSIDMw)

        for w in range(len(wids)-1):
            f.write(str(wids[w])+'\t')
        f.write(str(wids[-1])+'\n')

    f.close()
    plt.tight_layout()
    if withSIDM:
        plt.legend(['CDM','SIDM'])
        axs[0].savefig('/home/jw1624/H1-merian/figures/HIProfiles/HI_'+str(hID)+'_2.png')
    else:
        plt.savefig('/home/jw1624/H1-merian/figures/HIProfiles/HI_'+str(hID)+'.png')

# end makeHIprofile

# get haloIDs
cdmHalos, sidmHalos = util.getGalaxies()

print('Making HI Profiles')
# make figs
f=open('/home/jw1624/H1-merian/h1lines/widths.txt', 'w')
f.write('galaxy\tw50_cdm\tw20_cdm\tw10_cdm\tw50_sidm\tw20_sidm\tw10_sidm\n')
f.close()
for g in cdmHalos:
    print(' halo '+str(g)+'...', end='')
    if g in sidmHalos:
        makeHIprofile(g, withSIDM=True)
    makeHIprofile(g)
    print('done')
