import pynbody
import pylab
import numpy as np

from matplotlib import pyplot as plt

dpath = '/data/REPOSITORY/romulus_dwarf_zooms/r597.romulus25.3072g1HsbBH/'
f = 'r597.romulus25.3072g1HsbBH.004096'

s = pynbody.load(dpath+f)

# get first halo
h = s.halos()
h1 = h[1]

#print('ngas = %e, ndark = %e, nstar = %e\n'%(len(h1.gas),len(h1.dark),len(h1.star)))

# centering based on 1st halo
cen_pot = pynbody.analysis.halo.center(h1, mode='pot', retcen=True)
s['pos'] -= cen_pot

# a first image
s.physical_units()
pynbody.plot.image(h1.g, width=30, cmap='Blues', filename='figures/f1.png')

# side on
pynbody.analysis.angmom.sideon(h1, cen=(0,0,0))
pynbody.plot.image(h1.g, width = 30, cmap='Blues', filename='figures/f2.png')

# back to face
s.rotate_x(90)

# star profile
ps = pynbody.analysis.profile.Profile(h1.s, min=0.01, max=50, type='log')
pylab.clf()
pylab.plot(ps['rbins'], ps['density'])
pylab.seminology()
pylab.xlabel('$R$ [kpc]')
pylab.ylabel('$v_{circ}$ [km/s]')
pylab.legend()
