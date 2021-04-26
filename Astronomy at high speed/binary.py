#!/usr/bin/env python

"""Calculates eclipses of circular orbit binary, allowing for limb
darkening.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def format_axes(axes):
    """
    Fixes up ticks on plots into PGPLOT style
    """
    axes.tick_params(axis="x", direction="in")
    axes.tick_params(axis="y", direction="in")
    axes.tick_params(bottom=True, top=True, left=True, right=True)

def lcurve(phases, r1, r2, eps1, eps2, s1, s2, nr1, nr2, iangle):
    """Calculates the light curve at a set of orbital phases (defn: zero
    phase when star 2 is closest to Earth) of two stars in a circular
    orbit. It works by splitting the visible face of each star into a
    series of equal width annuli.

    Arguments:

       phases : numpy.ndarray
          array of phases. Phase 0 is defined as the time when star 2 is closest
          to Earth. 0-1 corresponds to one orbit.

       r1 : float
          radius of star 1, scaled by orbital separation = R1/a

       r2 : float
          radius of star 2, scaled by orbital separation = R2/a

       eps1 : float
          linear limb darkening coefficient of star 1. Set = 0 to
          turn off limb darkening. Range [0,1].

       eps2 : float
          linear limb darkening coefficient of star 2. Set = 0 to
          turn off limb darkening. Range [0,1].

       s1 : float
          surface brightness of star 1 seen at centre of visible face.

       s2 : float
          surface brightness of star 2 seen at centre of visible face.

       nr1 : int
          number of radii to represent star 1

       nr2 : int
          number of radii to represent star 2

       iangle : float
          orbital inclination, degrees

    """
    # a few checks on the inputs
    assert(r1 > 0 and r2 > 0)
    assert(eps1 >= 0 and eps1 <= 1 and eps2 >= 0 and eps2 <= 1)
    assert(s1 >= 0 and s2 >= 0)
    assert(nr1 > 0 and nr2 > 0)
    assert(iangle > 0 and iangle <= 90)

    cosi = np.cos(np.radians(iangle))
    sini = np.sin(np.radians(iangle))

    cphases = np.cos(2.*np.pi*phases)
    sphases = np.sin(2.*np.pi*phases)

    # array of scaled separations for each phase
    seps = np.sqrt(sphases**2+(cosi*cphases)**2)

    # radius and flux arrays. Visible faces of each star
    # are split into a series of equal width annuli. The central
    # radius and total flux from each annulus is saved.
    dr1 = r1/nr1
    r1s = np.arange(dr1/2,r1,dr1)
    fluxes1 = 2*np.pi*dr1*r1s*s1*(1-eps1+eps1*np.sqrt(1-(r1s/r1)**2))
    ftot1 = fluxes1.sum()

    dr2 = r2/nr2
    r2s = np.arange(dr2/2,r2,dr2)
    fluxes2 = 2*np.pi*dr2*r2s*s2*(1-eps2+eps2*np.sqrt(1-(r2s/r2)**2))
    ftot2 = fluxes2.sum()

    lc1 = np.empty_like(phases)
    lc2 = np.empty_like(phases)

    for n, (sep,cphase) in enumerate(zip(seps,cphases)):
        if cphase > 0:
            # star 2 is closest to Earth
            visibs = visible_fraction(r2, r1s, sep)
            lc1[n] = (fluxes1*visibs).sum()
            lc2[n] = ftot2
        else:
            # star 1 is closest to Earth
            visibs = visible_fraction(r1, r2s, sep)
            lc1[n] = ftot1
            lc2[n] = (fluxes2*visibs).sum()

    return (lc1+lc2,lc1,lc2)

def visible_fraction(rfront, rback, sep):
    """
    Computes visible fraction of circle(s).

    Arguments:

       rfront : float
          radius of obscuring circle

       rback : numpy.ndarray
          array of circles being obscured. Assumed concentric

       sep : float
          distance between centres of circles.
    """

    visib = np.empty_like(rback)

    # no obscuration
    on = sep >= rfront + rback
    visib[on] = 1.0

    # full obscuration
    off = rfront >= sep + rback
    visib[off] = 0.0

    # somewhere in between. compute cosines of angle subtended
    # relative to the line of centres in the background obscured
    # circles.
    mid = ~on & ~off
    cost = (rback[mid]**2+sep**2-rfront**2)/(2.*rback[mid]*sep)
    visib[mid] = 1-np.arccos(np.maximum(-1,np.minimum(1,cost)))/np.pi

    return visib


def solve_binary(phase3, phase4, f0, f1, f2):
    """Given values for the 3rd and 4th contact phases and the fluxes out
    of eclipse and at the middle of the primary and secondary eclipse,
    this solves for the scaled radii and surface brightness ratio of
    the two stars.

    Assumes: no limb darkening, total eclipses, star 2 > star 1 and
    star 1 eclipsed at phase 0.

    Arguments:

       phase3 : float
         phase of third contact. Phase 0 defined as star 2 eclipsing
         star 1. [units: cycles]

       phase4 : float
         phase of fourth contact [units: cycles]

       f0 : float
         flux out of eclipse

       f1 : float
         flux at mid-eclipse at phase 0

       f2 : float
         flux at mid-eclipse at phase 0.5

    Returns: (iangle,r1,r2,s1,s2)

    Orbital inclination in degrees, the scaled radii of each star and
    their surface brightnesses defined as the flux = visible area *
    surface brightness so that s1*Pi*r1**2 == flux from star 1 out of
    eclipse.

    See also: contacts_fluxes

    """

    assert(phase3 > 0 and phase4 > phase3)
    assert(phase4 < 0.25)
    assert(f0 > f1 and f0 > f2 and f2 > f1)

    # Compute radius ratio: alpha = r1/r2 < 1 by assumption.
    alpha = np.sqrt((f0-f2)/f1)

    # convert phases to radians
    phase3 *= 2.*np.pi
    phase4 *= 2.*np.pi

    # compute sines / cosines
    c3, s3 = np.cos(phase3), np.sin(phase3)
    c4, s4 = np.cos(phase4), np.sin(phase4)

    cosi = np.sqrt(
        (((1-alpha)*s4)**2-((1+alpha)*s3)**2) /
        (((1+alpha)*c3)**2-((1-alpha)*c4)**2)
    )
    r2 = np.sqrt(s4**2+(cosi*c4)**2)/(1+alpha)
    r1 = r2*alpha
    s2 = f1/(np.pi*r2**2)
    s1 = (f0-f1)/(np.pi*r1**2)
    iangle = np.degrees(np.arccos(cosi))

    return (iangle,r1,r2,s1,s2)

def contacts_fluxes(iangle,r1,r2,s1,s2):
    """Given the geometrical and surface brightness parameters of a
    binary, this returns values for the 3rd and 4th contact phases and
    the fluxes out of eclipse, and at the middle of the primary and
    secondary eclipse.

    Assumes: no limb darkening, total eclipses, star 2 > star 1, star 1
    has the highest surface brightness, and that star 1 is eclipsed at
    phase 0 (hence primary eclipse).

    Arguments:

       iangle : float
          orbital inclination [degrees, 0 < iangle <= 90]

       r1 : float
         radius of star 1 [units of orbital separation]

       r2 : float
         radius of star 2 [units of orbital separation]

       s1 : float
         surface brightness of star 1, defined so that the
         uneclipsed flux = Pi*r1**2*s1

       s2 : float
         surface brightness of star 2

    Returns: (phase3,phase4,f0,f1,f2)

    Third & fourth contact phases, and the fluxes out of eclipse and at
    mid primary and secondary eclipse.

    See also: solve_binary

    Notes :

    This can fail if invalid parameter combinations are supplied.

    """

    assert(iangle > 0 and iangle <= 90 and r2 >= r1 and s1 >= s2)
    assert(r1 > 0 and r2 >= r1 and s2 > 0 and s1 >= s2)
    cosi = np.cos(np.radians(iangle))
    assert(cosi < r2-r1)
    sini = np.sin(np.radians(iangle))

    sphi3 = np.sqrt((r2-r1)**2-cosi**2)/sini
    phase3 = np.arcsin(sphi3)/(2*np.pi)

    sphi4 = np.sqrt((r2+r1)**2-cosi**2)/sini
    phase4 = np.arcsin(sphi4)/(2*np.pi)

    f0 = np.pi*r1**2*s1 + np.pi*r2**2*s2
    f1 = np.pi*r2**2*s2
    f2 = np.pi*(r2**2-r1**2)*s2 + np.pi*r1**2*s1

    return (phase3, phase4, f0, f1, f2)

def func(x, phs, fs, fes):
    """
    Function to test a model vector of parameters against a set of data
    """
    iangle, r1, r2, s1, s2 = x
    nr1, nr2, eps1, eps2 = 200, 200, 0, 0
    try:
        lc, lc1, lc2 = lcurve(phs, r1, r2, eps1, eps2, s1, s2, nr1, nr2, iangle)
        return (((lc-fs)/fes)**2).sum()
    except:
        return float('Inf')

class Measure():
    """Cursor-based measurement of key features of an eclipsing binary
    light curve on a matplotlib plot. Uses "plt" as a synonym for
    matplotlib.pyplot.

    This prompts the user to define the five key observational
    measurements indicating the third & fourth contact phases along
    with the out of eclipse and mid-eclipse flux levels.

    At the end of the measurements the results (phase3, phase4, f0, f1, f2)
    are stored in an attribute of the Measure object called "results".
    """

    def __init__(self, ax):

        self.ax=ax
        self.results = []

        self.ndata = 0
        self.press=False
        self.move=False
        self.c1=self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmove)
        print('\n\nUse the mouse to define contact phases and flux levels.')
        print('You may want to pan and zoom the plot to help with your selections.')
        print('\n\nLeft-click to define the horizontal location of the third contact')

    def onclick(self,event):
        """
        Where stuff is done
        """
        if event.inaxes == self.ax and event.button == 1:
            if self.ndata == 0:
                self.results.append(event.xdata)
                print(f'.... stored third contact phase = {event.xdata}')
                print('\nLeft-click the location of the fourth contact')

            elif self.ndata == 1:
                self.results.append(event.xdata)
                print(f'.... stored fourth contact phase = {event.xdata}')
                print('\nLeft-click at the level of the out-of-eclipse flux')

            elif self.ndata == 2:
                self.results.append(event.ydata)
                print(f'.... stored out-of-eclipse flux = {event.ydata}')
                print('\nLeft-click at the level of the mid-primary-eclipse flux')

            elif self.ndata == 3:
                self.results.append(event.ydata)
                print(f'.... stored mid-primary-eclipse flux = {event.ydata}')
                print('\nLeft-click at the level of the mid-secondary-eclipse flux')

            elif self.ndata == 4:
                self.results.append(event.ydata)
                print(f'.... stored mid-secondary flux = {event.ydata}')
                plt.close()

            self.ndata += 1

    def onpress(self,event):
        self.press=True

    def onmove(self,event):
        if self.press:
            self.move=True

    def onrelease(self,event):
        # only call onclick in special circumstances
        if self.press and not self.move:
            self.onclick(event)
        self.press=False; self.move=False


class UpdateBinary:

    """
    Class to support use of matplotlib.animation.FuncAnimation
    for swift updating of plots to allow creation of movies.

    This one expects two axes one for shoing the binary, the
    other for the light curve.
    """
    
    def __init__(self, axb, axl, phases, lc, iangle, r1, r2, a1, a2):

        """
        Initialises the plot.
        """
        self.axb = axb
        self.axl = axl
        self.phases = phases
        self.lc = lc
        self.cosi = np.cos(np.radians(iangle))
        self.r1 = r1
        self.r2 = r2
        self.a1 = a1
        self.a2 = a2
        
        self.axb.set_xlim(-1.1,1.1)
        self.axb.set_ylim(-0.35,0.35)
        self.axb.set_aspect('equal')
        self.axb.set_xlabel('$x$')
        self.axb.set_ylabel('$y$')
        self.axb.set_title(f'i = {iangle}, r1 = {r1}, r2 = {r2}')
        
        self.axl.set_xlim(-0.1,0.9)
        self.axl.set_ylim(0,0.6)
        self.axl.set_xlabel('Orbital phase [cycles]')
        self.axl.set_ylabel('Flux')

        # Create and plot stars
        cosp = np.cos(2.*np.pi*self.phases[0])
        sinp = np.sin(2.*np.pi*self.phases[0])
        self.star1 = Circle((-self.a1*sinp, self.a1*self.cosi*cosp), self.r1, color='b', zorder=5)
        if cosp > 0:
            self.star2 = Circle((self.a2*sinp, -self.a2*self.cosi*cosp), self.r2, color='r', zorder=10, alpha=0.9)
        else:
            self.star2 = Circle((self.a2*sinp, -self.a2*self.cosi*cosp), self.r2, color='r', zorder=0)

        self.axb.add_artist(self.star1)
        self.axb.add_artist(self.star2)

        # Plot light curve
        self.axl.plot(self.phases, self.lc, '0.7', zorder=0)
        self.lcpl, = self.axl.plot(self.phases[:1], self.lc[:1], 'g', lw=3, zorder=1)
        self.lcpt, = self.axl.plot(self.phases[0], self.lc[0], 'ok')
        

    def __call__(self, n):
        
        # Update stars
        cosp = np.cos(2.*np.pi*self.phases[n])
        sinp = np.sin(2.*np.pi*self.phases[n])
        self.star1.set_center((-self.a1*sinp, self.a1*self.cosi*cosp))
        self.star2.set_center((self.a2*sinp, -self.a2*self.cosi*cosp))
        if cosp > 0:
            self.star2.set_zorder(10)
        else:
            self.star2.set_zorder(0)

        # Update lc plot
        self.lcpl.set_data(self.phases[:n], self.lc[:n])
        self.lcpt.set_data(self.phases[n], self.lc[n])

        # Return list of "artist" objects that have been updated
        # This is the essential property required of the "func"
        # argument of FuncAnimation
        return (self.star1, self.star2, self.lcpl, self.lcpt)        
                
if __name__ == '__main__':

    phase3, phase4, f0, f1, f2 = 0.0073970705738562945, 0.0664016960580229, 1.4398232550213454, 0.2304931639141266, 1.3277282161981918

    print(solve_binary(phase3, phase4, f0, f1, f2))
