# util_galaxies.py
# Joseph Wick
# H1-merian
#
# utility functions for galaxy analysis
#

import numpy as np
import random

class util_galaxies:

    # compute_CoM
    # computes center of mass given arrays of positions and masses
    @staticmethod
    def compute_CoM(pos, mass):
        '''
        compute_CoM()

        Computes center of mass given positions and masses of particles. Used
        when there is not a halo catalog for a given sim file, otherwise it is
        better to use Pynbody centering (momentum)

        Parameters
        ----------
        pos : array_like
            positions of particles

        mass : array_like
            masses of particles

        Returns
        -------
        cen : array_like
            Center of mass

        '''
        pos_np = np.array(pos)
        mass_np = np.array(mass)

        mtot = mass_np.sum()
        cen = np.sum(mass_np * pos_np.transpose(), axis=1) / mtot

        return cen

    # compute_vdisp_wtd
    # compute mass weighted velocity dispersion
    @staticmethod
    def compute_vdisp_global(vel_all, mass_all, vel_t, weights=1):
        '''
        compute_vdisp_global

        Computes weighted velocity dispersion using standard deviation.
        Argument `weights` defaults to 1 for unweighted functionality.

        Parameters
        -----
        vel_all : array_like
            velocities of all stars. Used for center of mass velocity

        mass_all: array_like
            masses of all stars. Used for center of mass velocity

        vel_t : array_like
            velocities of target stars we want to calculate v disp for

        weights : array_like
            weights of target stars we want to calculate vdisp for

        Returns
        -----
        v_disp : float
            mass weighted standard deviation of difference from systemic velocity
            for stars in vel_t and weights
        '''

        # check that both `_all` arrays have the same length
        if len(vel_all)!=len(mass_all):
            print('Error: lengths of `vel_all` and `mass_all` must be identical')
            sys.exit()
        # check that both `_t` arrays have the same length
        if len(vel_all)!=len(mass_all):
            print('Error: lengths of `vel_t` and `weights` must be identical')
            sys.exit()

        # center of mass velocity / systemic velocity
        v_CoM = np.sum(vel_all * mass_all[:, None], axis=0) / np.sum(mass_all)

        # difference from systemic velocity
        v_diffs = vel_t - v_CoM
        # take magnitude
        v_diffs_tot = np.linalg.norm(v_diffs, axis=1)

        # mass weighted difference from mean vdiff
        vdisp = np.sum( ((v_diffs_tot-np.mean(v_diffs_tot))**2) * weights)
        # divide by total of weights
        vdisp = vdisp/np.sum(weights)
        # sqrt
        vdisp = np.sqrt(vdisp)

        return float(vdisp)

    # compute_vdisp_los()
    # computes median line of sight velocity dispersion
    def compute_vdisp_los(vel_all, mass_all, pos_t, vel_t, hmr, weights=[]):
    '''
    compute_vdisp_los()

    Computes (weighted) line of sight velocity dispersion. Weights default to one
    for nonweighted functionality if desired.

    Calculates line of sight velocity dispersion at 100 random observation points
    and returns the median.

    Parameters
    ----------
    vel_all : array_like
        Velocities for all stars. Used for systematic velocity subtraction.

    mass_all : array_like
        Masses of all stars. Used for systematic velocity subtraction.

    pos_t : array_like
        Centered positions of stars that vdisp are calculated for.

    vel_t : array_like
        velocities of stars that vdisp are calculated for.

    hmr : float
        Half mass radius of selected galaxiy.

    weights : array_like (optional)
        Optional weights for vdisp calculation. Defaults to all ones.

    Returns
    -------

    vdisp_los : float
        Median (weighted) velocity dispersion from 100 random line of sight
        calculations.
    '''

    R = 10*hmr

    if len(weights == 0):
        weights = np.ones(len(pos_t))

    # subtract systematic/CoM velocity
    v_CoM = np.sum(vel_all * mass_all[:, None], axis=0) / np.sum(mass_all)
    vdiffs = np.array(vel_t - v_CoM)
    #vdiffs_mag = np.linalg.norm(vdiffs, axis=1)

    allvdisps = []
    for i in range(100):
        # select viewing location
        theta = random.uniform(0, np.pi)
        phi = random.uniform(0, 2*np.pi)

        x = R*np.sin(theta)*np.cos(phi)
        y = R*np.sin(theta)*np.sin(phi)
        z = R*np.cos(theta)

        # unit vector from particle position to observation point
        a = np.array([x,y,z]-pos_t)
        a = a/np.linalg.norm(a, axis=1)[:,None]
        # projection of velocity along vector a
        vel_los = np.linalg.norm(vdiffs*a, axis=1)

        weightedstd = np.sum( ((vel_los-np.mean(vel_los))**2) * weights)
        weightedstd = weightedstd/np.sum(weights)
        weightedstd = np.sqrt(weightedstd)

        allvdisps.append(weightedstd)

    return np.median(allvdisps)

    # compute_massRadius()
    # computes radius within percentage of mass lie
    def compute_massRadius(positions, masses, outerR, acc, frac=0.5, maxiter=10000):
        '''
        compute_massRadius

        Computes radius that encloses frac*sum(masses) amount of matter.
        Defaults to half mass radius.

        Parameters
        ----------

        positions : array_like
            Particle positions

        masses : array_like
            Particle masses

        outerR : float
            Maximum radius to check

        acc : float
            Percent accuracy expressed as a decimal (ie 0.01 = 1% error tolerance)

        frac : float (optional)
            Fraction of mass to enclose, expressed as a decimal (0.5 = 50%).
            Defaults to 0.5

        maxiter : float (optional)
            Maximum bisections to perform. Defaults to 1000000.

        Returns
        -------
        r : float
            Fractional mass radius. Returns -1 if maxiter reached.
        '''

        innerLim = 0.0
        outerLim = outerR

        # center by mass
        mtot = masses.sum()
        cen = np.sum(masses * positions.transpose(),
             axis=1) / mtot
        cenpos = positions - cen

        # take radii
        pRadii = np.linalg.norm(cenpos, ord=2, axis=1)
        mTot = np.sum(masses)

        r = outerR/2
        hm = np.sum(masses[pRadii < r])
        n=0
        while(hm < (frac-acc)*mTot or hm > (frac+acc)*mTot):
            if hm > 0.5*mTot: # too big, decrease r
                outerLim = r
                r = innerLim + (outerLim-innerLim)/2
            elif hm < 0.5*mTot: # too small, increase r
                innerLim = r
                r = innerLim + (outerLim-innerLim)/2
            hm = sum(masses[pRadii < r])

            n += 1
            if n>maxiter:
                return -1
                break

        return r
    #
