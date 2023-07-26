# makeGalImages.py
# also exports median v_disp for each galaxy

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import pynbody
import numpy as np
import pynbody.plot.sph as sph
import matplotlib.pylab as plt


# get current halos
cdmHalos, sidmHalos, adiabaticHalos = util.getGalaxies()

# image function
def makeStarImage(hID, width=20, sidm=0):
    fileadd = '.romulus25.3072g1HsbBH.004096'
    if sidm:
        fileadd = '.romulus25cvdXsec.3072g1HsbBH.004096'

    v_disp = -1

    fig = plt.figure(figsize=(8,8), facecolor='w')

    f = util.getfilepath(hID)[sidm]+'/r'+str(hID)+fileadd

    s = pynbody.load(f)
    s.physical_units()

    h = s.halos()
    v_disp = np.median(np.array(h[1].g['v_disp']))

    # center on the largest halo and align the disk
    pynbody.analysis.angmom.faceon(h[1])
    pynbody.plot.stars.render(s, width=str(width)+' kpc')

    tag = 'CDM'
    if sidm: tag = 'SIDM'
    plt.title('r'+str(hID) + ' '+tag, fontsize=20, fontfamily='serif')

    plt.tight_layout()

    if sidm:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(hID)+'SIDMStarImg_f.png'
        plt.savefig(fname)
    else:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(hID)+'CDMStarImg_f.png'
        plt.savefig(fname)

    # side on
    fig = plt.figure(figsize=(8,8), facecolor='w')
    pynbody.analysis.angmom.sideon(h[1])
    pynbody.plot.stars.render(s, width=str(width)+' kpc')

    plt.title('r'+str(hID) + ' '+tag, fontsize=20, fontfamily='serif')

    plt.tight_layout()

    if sidm:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(hID)+'SIDMStarImg_s.png'
        plt.savefig(fname)
    else:
        fname = '/home/jw1624/H1-merian/figures/images/r'+str(hID)+'CDMStarImg_s.png'
        plt.savefig(fname)

    # write v_disp
    if sidm:
        f = open('/home/jw1624/H1-merian/csvs/vDisp_sidm.txt', 'a')
        f.write(str(hID)+','+str(v_disp)+'\n')
        f.close()
    else:
        f = open('/home/jw1624/H1-merian/csvs/vDisp_cdm.txt', 'a')
        f.write(str(hID)+','+str(v_disp)+'\n')
        f.close()

    print('Finished halo '+str(hID) + ' in '+tag)
# end make StarImage

# now actually make the image
f = open('/home/jw1624/H1-merian/csvs/vDisp_cdm.txt', 'w')
f.write('galaxy,v_disp\n')
f.close()
f = open('/home/jw1624/H1-merian/csvs/vDisp_sidm.txt', 'w')
f.write('galaxy,v_disp\n')
f.close()
for cdmg in cdmHalos:
    makeStarImage(cdmg)
#for sidmg in sidmHalos: # only want CDM currently
#    makeStarImage(sidmg, sidm=1)
