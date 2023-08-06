from util_os import util_os

import pynbody



##

# test retrieval of z=0
cdmGals = util_os.getGalaxies()[0]
print('Testing z=0 retrieval...', end='')
for gal in cdmGals:
    z0 = util_os.getfilepath_cdm(gal, 0)
    sCDM = pynbody.load(z0)
print('Success\n')

# test retrieval of all CDM timesteps
for gal in cdmGals:
    print('Testing all retrieval for r'+str(gal)+'...', end='')
    nts = util_os.getNumTimesteps(gal)
    for ts in range(nts):
        f = util_os.getfilepath_cdm(gal, ts)
    print('Success')
    print(nts)
