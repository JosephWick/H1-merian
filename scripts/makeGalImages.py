import sys
sys.path.insert(0, '../util/')
from util import util

import pynbody
import pynbody.plot.sph as sph
import matplotlib.pylab as plt


# get current halos 
cdmHalos, sidmHalos = getGalaxies()

# image function 
def makeStarImage(hID, width=10, sidm=0):
    fileadd = '.romulus25.3072g1HsbBH.004096'
    if sidm:
        fileadd = '.romulus25cvdXsec.3072g1HsbBH.004096'

    fig = plt.figure(figsize=(8,8), facecolor='w')

    f = getfilepath(hID)[sidm]+'/r'+str(hID)+fileadd

    s = pynbody.load(f)
    s.physical_units()

    h = s.halos()

    # center on the largest halo and align the disk
    pynbody.analysis.angmom.faceon(h[1])

    pynbody.plot.stars.render(s,width=str(width)+' kpc')

    tag = 'CDM'
    if sidm: tag = 'SIDM'
    plt.title('r'+str(hID) + ' '+tag, fontsize=20, fontfamily='serif')

    if sidm:
        plt.savefig('../figures/images/r'+str(hID)+'SIDMStarImg.png')
    else:
        plt.savefig('../figures/images/r'+str(hID)+'CDMStarImg.png')
# end make StarImage

# now actually make the image 
for cdmg in cdmHalos:
    makeStarImage(cdmg)
for sidmg in sidmHalos:
    makeStarImage(sidmg)


