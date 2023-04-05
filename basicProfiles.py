import pynbody
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import pyplot as plt

# note these are all CDM
files = ['/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r431.romulus25.3072g1HsbBH/r431.romulus25.3072g1HsbBH.004096/r431.romulus25.3072g1HsbBH.004096',
         '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r468.romulus25.3072g1HsbBH/r468.romulus25.3072g1HsbBH.004096/r468.romulus25.3072g1HsbBH.004096',
         '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r492.romulus25.3072g1HsbBH/r492.romulus25.3072g1HsbBH.004096/r492.romulus25.3072g1HsbBH.004096',
         '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r502.romulus25.3072g1HsbBH/r502.romulus25.3072g1HsbBH.004096/r502.romulus25.3072g1HsbBH.004096',
         '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r569.romulus25.3072g1HsbBH/r569.romulus25.3072g1HsbBH.004096/r569.romulus25.3072g1HsbBH.004096']
haloIDs = [431, 468, 492, 502, 569]

for i,f in enumerate(files):
    print('Working on halo '+str(haloIDs[i]))

    s = pynbody.load(f)

    # get first halo
    h = s.halos()
    h1 = h[1]

    #print('ngas = %e, ndark = %e, nstar = %e\n'%(len(h1.gas),len(h1.dark),len(h1.star)))

    # centering based on 1st halo
    cen_pot = pynbody.analysis.halo.center(h1, mode='pot', retcen=True)
    s['pos'] -= cen_pot

    # a first image
    s.physical_units()
    pynbody.plot.image(h1.g, width=30, cmap='Blues', filename='figures/r'+str(haloIDs[i])+'CDM_img1.png')

    # side on
    pynbody.analysis.angmom.sideon(h1, cen=(0,0,0))
    pynbody.plot.image(h1.g, width = 30, cmap='Blues', filename='figures/r'+str(haloIDs[i])+'CDM_img2.png')

    # back to face
    s.rotate_x(90)

    # profile range
    pmin = 0.01
    #pmax = pynbody.analysis.halo.virial_radius(h1) / 4
    pmax = 4

    # rotation curve
    plt.figure()
    pd = pynbody.analysis.profile.Profile(h1.d, rmin=pmin, rmax=pmax, type='lin')
    pg = pynbody.analysis.profile.Profile(h1.g, rmin=pmin, rmax=pmax, type='lin')
    ps = pynbody.analysis.profile.Profile(h1.s, rmin=pmin, rmax=pmax, type='lin')
    p  = pynbody.analysis.profile.Profile(h1,   rmin=pmin, rmax=pmax, type='lin')

    for prof, name in zip([p,pd,ps,pg], ['total', 'dm', 'stars', 'gas']):
        plt.plot(prof['rbins'], prof['v_circ'], label=name)

    plt.title('Rotations')
    plt.xlabel('$R$ [kpc]')
    plt.ylabel('$v_{circ}$ [km/s]')

    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.tight_layout()
    plt.savefig('figures/r'+str(haloIDs[i])+'_CDM_rotation.png')

    # density profiles
    pmax = pynbody.analysis.halo.virial_radius(h1) / 4

    plt.figure()
    for prof, name in zip([p,pd,ps,pg], ['total', 'dm', 'stars', 'gas']):
        plt.plot(prof['rbins'], prof['density'], label=name)

    plt.title('Density')
    plt.xlabel('$R$ [kpc]')
    plt.ylabel('density [M$_{\odot}$ kpc$^{-3}$]')

    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.tight_layout()
    plt.savefig('figures/r'+str(haloIDs[i])+'_CDM_density.png')
