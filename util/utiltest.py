from util_os import util_os

import pynbody



##

cdmGals = util_os.getGalaxies()[0]
print(cdmGals)
for gal in currentGals:
    z0 = util_os.getfilepath(gal,0)
    sCDM = pynbody.load(z0)
