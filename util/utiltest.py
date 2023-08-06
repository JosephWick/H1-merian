from util_os import util_os

import pynbody



##

# test retrieval of z=0
cdmGals = util_os.getGalaxies()[0]
print('Testing z=0 retrieval...')
for gal in cdmGals:
    z0 = util_os.getfilepath_cdm(gal, 0)
    sCDM = pynbody.load(z0)
print('success')
