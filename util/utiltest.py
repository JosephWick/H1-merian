from util_os import util_os

immport pynbody



##

currentGals = util_os.getGalaxies()
for gal in currentGals:
    z0 = util_os.getfilepath(gal,0)
    sCDM = pynbody.load(zo)
