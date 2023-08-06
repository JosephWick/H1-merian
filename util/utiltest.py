from util_os import util_os

import pynbody



##

currentGals = util_os.getGalaxies()
for gal in currentGals:
    z0 = util_os.getfilepath(gal,0)
    sCDM = pynbody.load(z0)
