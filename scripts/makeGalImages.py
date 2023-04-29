import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import pynbody
import pynbody.plot.sph as sph
import matplotlib.pylab as plt


# get current halos
cdmHalos, sidmHalos, adiabaticHalos = util.getGalaxies()

# image function
def makeStarImage(hID, width=15, sidm=0):
    fileadd = '.romulus25.3072g1HsbBH.004096'
    if sidm:
        fileadd = '.romulus25cvdXsec.3072g1HsbBH.004096'

    fig, axs = plt.subplots(2,2, figsize=(8,8), facecolor='w')

    f = getfilepath(hID)[sidm]+'/r'+str(hID)+fileadd

    s = pynbody.load(f)
    s.physical_units()

    h = s.halos()

    # center on the largest halo and align the disk
    pynbody.analysis.angmom.faceon(h[1])
    pynbody.plot.stars.render(s, subplot=axs[0], width=str(width)+' kpc')

    # side on
    pynbody.analysis.angmom.sideon(h[1])
    pynbody.plot.stars.render(s, subplot=axs[1], width=str(width)+' kpc')

    tag = 'CDM'
    if sidm: tag = 'SIDM'
    plt.title('r'+str(hID) + ' '+tag, fontsize=20, fontfamily='serif')

    plt.tight_layout()

    if sidm:
        plt.savefig('/home/jw1624/H1-merian/figures/images/r'+str(hID)+'SIDMStarImg.png')
    else:
        plt.savefig('/home/jw1624/H1-merian/figures/images/r'+str(hID)+'CDMStarImg.png')
# end make StarImage

# now actually make the image
for cdmg in cdmHalos:
    makeStarImage(cdmg)
for sidmg in sidmHalos:
    makeStarImage(sidmg)
