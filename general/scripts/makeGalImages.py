# makeGalImages.py
# also exports median v_disp for each galaxy

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util_os import util_os

import pynbody
import numpy as np
import pynbody.plot.sph as sph
import matplotlib.pylab as plt


# get current halos
cdmHalos, sidmHalos, adiabaticHalos = util_os.getGalaxies()

# image function
def makeStarImage(gal, width=20, sidm=0):
    f = util_os.getfilepath_cdm(gal, 0)
    if sidm:
        f = util_os.getfilepath_sidm(gal, 0)

    v_disp = -1

    fig = plt.figure(figsize=(8,8), facecolor='w')

    s = pynbody.load(f)
    s.physical_units()

    h = s.halos()
    v_disp = np.median(np.array(h[1].g['v_disp']))

    # center on the largest halo and align the disk
    pynbody.analysis.angmom.faceon(h[1])
    pynbody.plot.stars.render(s, width=str(width)+' kpc')

    tag = 'CDM'
    if sidm: tag = 'SIDM'
    plt.title('r'+str(gal) + ' '+tag, fontsize=20, fontfamily='serif')

    plt.tight_layout()

    if sidm:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(gal)+'SIDMStarImg_f.png'
        plt.savefig(fname)
    else:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(gal)+'CDMStarImg_f.png'
        plt.savefig(fname)

    # side on
    fig = plt.figure(figsize=(8,8), facecolor='w')
    pynbody.analysis.angmom.sideon(h[1])
    pynbody.plot.stars.render(s, width=str(width)+' kpc',  resolution=1000)

    plt.title('r'+str(gal) + ' '+tag, fontsize=20, fontfamily='serif')

    plt.tight_layout()

    if sidm:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(gal)+'SIDMStarImg_s.png'
        plt.savefig(fname, dpi=600)
    else:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(gal)+'CDMStarImg_s.png'
        plt.savefig(fname, dpi=600)

    # write v_disp
    if sidm:
        f = open('/home/jw1624/H1-merian/csvs/vDisp_sidm.txt', 'a')
        f.write(str(gal)+','+str(v_disp)+'\n')
        f.close()
    else:
        f = open('/home/jw1624/H1-merian/csvs/vDisp_cdm.txt', 'a')
        f.write(str(gal)+','+str(v_disp)+'\n')
        f.close()

    print('Finished halo '+str(gal) + ' in '+tag)
    plt.close()
# end make StarImage

# sizes
widths = [20, 20, 10, 20, 10, 20, 20, 10, 15, 20, 20, 15, 15, 15, 10, 20]

# now actually make the image
f = open('/home/jw1624/H1-merian/csvs/vDisp_cdm.txt', 'w')
f.write('galaxy,v_disp\n')
f.close()
f = open('/home/jw1624/H1-merian/csvs/vDisp_sidm.txt', 'w')
f.write('galaxy,v_disp\n')
f.close()
for i,cdmg in enumerate(cdmHalos):
    makeStarImage(cdmg, width=widths[i])
#for sidmg in sidmHalos: # only want CDM currently
#    makeStarImage(sidmg, sidm=1)