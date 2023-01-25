import numpy as np
from glob import glob
from astropy.io import fits
import matplotlib.pyplot as plt

# data
datadir = '/data/REPOSITORY/romulus_dwarf_zooms/r597.romulus25.3072g1HsbBH'
datafiles = glob(datadir+'/r597.alfalfa_hi.*.fits')

print(datafiles)

for fname in datafiles:
    data = fits.open(fname)[0].data
    ftitle = fname[-17:-4]

    # define square to integrate over
    start = 0
    end = len(data[0])

    # --- INTEGRATE ---

    # integrated values array
    I = np.nan_to_num(np.array(data[:,start,start]), True, nan=0.0)

    for i in range(start+1,end):
        d = np.array(data[:,i,i])
        d = np.nan_to_num(d, False,  nan=0.0)
        I += d

    vmax = max(I)

    # --- MAKE FIGURE ---
    plt.figure(figsize=(12,8))
    plt.plot(I, linewidth=4)

    # prettify
    plt.xlabel('velocity [km/s]')
    plt.ylabel('intensity [Jy?]')
    plt.title(fname[len(datadir):])

    # add in percentages of max
    intervals = [10, 20, 50]
    c = ['red', 'orange', 'yellow']
    for i,p in enumerate(intervals):
        val = (p/100)*vmax
        plt.plot([0,35],[val,val], linewidth=1, c=c[i])

    plt.legend(['intensity', 'W10', 'W20', 'W50'])

    plt.savefig('figures/prelimH1_'+ftitle)
