from util_os import util_os

import pynbody



##

cdmGals = util_os.getGalaxies()[0]
for gal in cdmGals:
    print('Retrieving z=0 r'+str(gal)+'...', end='')
    z0 = util_os.getfilepath(gal,0)
    sCDM = pynbody.load(z0)
    print('done')
