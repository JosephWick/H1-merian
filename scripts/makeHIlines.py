import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import glob


# mu_n for kurtosis (defn from El Bhadri et al 2018)
def mu_n(v,s, n):
    # s is line profile
    # v is velocity bins
    vbar = np.array(s*v).mean()
    print('vbar: '+str(vbar))

    top = (s*(v-vbar)**n).sum()
    bottom = s.sum()
    return top/bottom

# getKurtosis
def getKurtosis(v,s):
    return mu_n(v,s,4)/mu_n(v,s,3) - 3

def makeHIprofile(hID, withSIDM=False, doExport=True):
    f=open('/home/jw1624/H1-merian/csvs/HI_widths.txt', 'a')

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
    K_cdm = -1
    for i in range(len(fcdm)):
        cdmx = pd.read_csv(fcdm[i], sep='\s+', header=None)[0]
        cdmy = pd.read_csv(fcdm[i], sep='\s+', header=None)[1]

        K_cdm = getKurtosis(np.array(cdmx), np.array(cdmy))

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

            axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cCDMw, label='_nolegend_')
            axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2, marker=markers[j], c=cCDMw)

        # do line widths of sidm
        K_sidm = -1
        if withSIDM:
            sidmx = pd.read_csv(fsidm[i], sep='\s+', header=None)[0]
            sidmy = pd.read_csv(fsidm[i], sep='\s+', header=None)[1]

            axs[i].plot(sidmx, sidmy, linewidth=lw-1, c=cSIDM)

            K_sidm = getKurtosis(np.array(sidmx), np.array(sidmy))

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

                axs[i].plot([x1,x2],[val,val], linewidth=lwW, zorder=2, c=cSIDMw, label='_nolegend_')
                axs[i].scatter([x1,x2],[val,val], s=25, label='_nolegend_', zorder=2, marker=markers[j], c=cSIDMw)

        if doExport:
            f.write(str(hID)+'_'+orientations[i]+',')
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

            # kurtosis
            f.write(str(K_cdm)+','+str(K_sidm)+'\n')

            # mean v_disp



    f.close()
    plt.tight_layout()
    if withSIDM:
        axs[0].legend(['CDM','SIDM'])
        plt.savefig('/home/jw1624/H1-merian/figures/HIProfiles/HI_'+str(hID)+'_2.png')
    else:
        plt.savefig('/home/jw1624/H1-merian/figures/HIProfiles/HI_'+str(hID)+'.png')

# end makeHIprofile

# get haloIDs
cdmHalos, sidmHalos, adiabaticHalos = util.getGalaxies()

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
    if g in sidmHalos:
        makeHIprofile(g, withSIDM=True, doExport=True)
        makeHIprofile(g, withSIDM=False,doExport=False)
    else:
        makeHIprofile(g, withSIDM=False,doExport=True)
    print('done')
    break
