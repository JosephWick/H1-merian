import '../pynbody'
import pylab


dpath = '/data/REPOSITORY/romulus_dwarf_zooms/r597.romulus25.3072g1HsbBH/'
f = 'r597.romulus25.3072g1HsbBH.004096'

s = pynbody.load(dpath+f)
print('loaded')
