#!/usr/bin/env python

"""
Plots the ellipse geometry of star 2 relative to star 1 in a couple of ways.

Plots for the notebook text.
"""

import numpy as np
import matplotlib.pyplot as plt
import binary

if __name__ == '__main__':

    # Define binary
    iangle = 87
    r1 = 0.15
    r2 = 0.25
    s1 = 5.0
    s2 = 1.5

    # Particular phase to place star 2
    phase = 0.05

    cosi = np.cos(np.radians(iangle))
    assert(cosi < r2-r1)

    # Compute path of star 2 splitting into front and back
    # sections to give a bit of "3D ness"
    front = np.linspace(-np.pi/2,np.pi/2,500)
    back = np.linspace(np.pi/2,3*np.pi/2,500)
    xf =  np.sin(front)
    yf = -cosi*np.cos(front)
    xb =  np.sin(back)
    yb = -cosi*np.cos(back)

    # First plot to illustrate ellipse
    fig, ax = plt.subplots(figsize=(10,5))
    binary.format_axes(ax)
    ax.plot(xf,yf,'r--',zorder=10)
    ax.plot(xb,yb,'r--',zorder=0)
    ax.plot([-1.1,1.1],[0,0],'--',color='0.7',zorder=5, alpha=0.5)
    ax.plot([0,0],[-0.3,0.3],'--',color='0.7',zorder=5, alpha=0.5)

    star1 = plt.Circle((0., 0.), r1, color='b', zorder=3)

    phi = 2.*np.pi*phase
    x2,y2 = np.sin(phi), -cosi*np.cos(phi)
    star2 = plt.Circle((x2,y2), r2, color='r', zorder=15)

    ax.plot(0,0,'k.',zorder=20)
    ax.plot(x2,y2,'k.',zorder=20)
    ax.add_artist(star1)
    ax.add_artist(star2)
    ax.set_aspect('equal')
    ax.set_ylim(-0.31,0.31)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    plt.title(f'$i = {iangle}^\circ$, $r_1 = {r1}$, $r_2 = {r2}$, $\phi = {phase}$')
    plt.savefig('eclipse-geometry1.png',bbox_inches='tight')
    plt.close()

    # Second plot to show contact phases
    fig, axs = plt.subplots(4,2,figsize=(10,6.8),sharex=True,sharey=True)

    phase3,phase4,f0,f1,f2 = binary.contacts_fluxes(iangle,r1,r2,s1,s2)

    phases = [-phase4,-phase3, phase3, phase4]

    for n, (phase, ax) in enumerate(zip(phases, axs[:,0])):

        binary.format_axes(ax)
        ax.plot(xf,yf,'--',color='0.5',zorder=10, alpha=0.5)
        ax.plot(xb,yb,'--',color='0.5',zorder=0, alpha=0.5)

        star1 = plt.Circle((0., 0.), r1, color='b', zorder=3)
        phi = 2.*np.pi*phase
        x2,y2 = np.sin(phi), -cosi*np.cos(phi)
        star2 = plt.Circle((x2,y2), r2, color='r', zorder=15, alpha=0.9)

        ax.add_artist(star1)
        ax.add_artist(star2)
        ax.set_aspect('equal')
        ax.set_ylim(-0.31,0.31)
        ax.text(1.0,0.2,f'{n+1}')
        ax.set_ylabel('y')
        if n == len(axs)-1:
            ax.set_xlabel('$x$')

    phases = [0.5-phase4,0.5-phase3, 0.5+phase3, 0.5+phase4]

    for n, (phase, ax) in enumerate(zip(phases, axs[:,1])):

        binary.format_axes(ax)
        ax.plot(xf,yf,'--',color='0.5',zorder=10, alpha=0.5)
        ax.plot(xb,yb,'--',color='0.5',zorder=0, alpha=0.5)

        star1 = plt.Circle((0., 0.), r1, color='b', zorder=3)
        phi = 2.*np.pi*phase
        x2,y2 = np.sin(phi), -cosi*np.cos(phi)
        star2 = plt.Circle((x2,y2), r2, color='r', zorder=1)

        ax.add_artist(star1)
        ax.add_artist(star2)
        ax.set_aspect('equal')
        ax.set_ylim(-0.35,0.35)

        ax.text(1.0,0.2,f'{n+1}')

        if n == len(axs)-1:
            ax.set_xlabel('x')

    plt.tight_layout()
    plt.savefig('eclipse-geometry2.png',bbox_inches='tight')
    plt.close()

    # Third plot to show the light curve
    fig = plt.figure(figsize=(10,7))
    grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.2)

    ax1 = fig.add_subplot(grid[0, :])
    binary.format_axes(ax1)
    ax2 = fig.add_subplot(grid[1, 0], sharey=ax1)
    binary.format_axes(ax2)
    ax3 = fig.add_subplot(grid[1, 1], sharey=ax1)
    binary.format_axes(ax3)

    p1, p2, nph = -0.15, 1.65, 2000
    phases = np.linspace(p1, p2, nph)
    lc,lc1,lc2 = binary.lcurve(phases, r1, r2, 0, 0, s1, s2, 1000, 1000, iangle)

    ax1.plot(phases, lc, 'b')
    ax1.set_ylim(0,0.8)
    ax1.set_ylabel('Flux')

    # Arrows to indicate key flux values
    ax1.arrow(
        0,0,0,f1,
        width=0.015, head_width=0.04, head_length=0.1,
        length_includes_head=True,fc='g',ec=None,lw=0,
        alpha=0.5
    )
    ax1.text(0.02,0.08,'$f_1$',fontsize=16)

    ax1.arrow(
        0.25,0,0,f0,
        width=0.015, head_width=0.04, head_length=0.1,
        length_includes_head=True,fc='g',ec=None,lw=0,
        alpha=0.5
    )
    ax1.text(0.272,0.08,'$f_0$',fontsize=16)

    ax1.arrow(
        0.5,0,0,f2,
        width=0.015, head_width=0.04, head_length=0.1,
        length_includes_head=True,fc='g',ec=None,lw=0,
        alpha=0.5
    )
    ax1.text(0.522,0.08,'$f_2$',fontsize=16)
    ax1.set_xlabel('Orbital phase [cycles]')
    ax1.set_xlim(p1,p2)

    llc,llc1,llc2 = binary.lcurve(phases, r1, r2, 0.6, 0.6, s1, s2, 1000, 1000, iangle)
    llc *= lc.max()/llc.max()

    # zoom on the primary eclipse
    width, off, frac = 0.09, 0.002, 0.15
    y1,ymax = ax2.set_ylim()

    ax2.plot([-phase4,-phase4],[0,ymax],'--',color='0.7')
    ax2.text(-phase4+off, frac*ymax, '$\phi_1$',fontsize=16)
    ax2.plot([-phase3,-phase3],[0,ymax],'--',color='0.7')
    ax2.text(-phase3+off, frac*ymax, '$\phi_2$',fontsize=16)
    ax2.plot([phase3,phase3],[0,ymax],'--',color='0.7')
    ax2.text(phase3+off, frac*ymax, '$\phi_3$',fontsize=16)
    ax2.plot([phase4,phase4],[0,ymax],'--',color='0.7')
    ax2.text(phase4+off, frac*ymax, '$\phi_4$',fontsize=16)

    ax2.plot(phases, llc, '0.7')
    ax2.plot(phases, lc, 'b')
    ax2.set_xlim(-width,width)
    ax2.set_ylabel('Flux')
    ax2.set_xlabel('Orbital phase [cycles]')

    # zoom on the secondary eclipse
    ax3.plot(phases, llc, '0.7')
    ax3.plot(phases, lc, 'b')
    ax3.set_xlim(0.5-width, 0.5+width)
    ax3.set_xlabel('Orbital phase [cycles]')

    ax3.plot([0.5-phase4,0.5-phase4],[0,ymax],'--',color='0.7')
    ax3.text(0.5-phase4+off, frac*ymax, "$\phi'_1$",fontsize=16)
    ax3.plot([0.5-phase3,0.5-phase3],[0,ymax],'--',color='0.7')
    ax3.text(0.5-phase3+off, frac*ymax, "$\phi'_2$",fontsize=16)
    ax3.plot([0.5+phase3,0.5+phase3],[0,ymax],'--',color='0.7')
    ax3.text(0.5+phase3+off, frac*ymax, "$\phi'_3$",fontsize=16)
    ax3.plot([0.5+phase4,0.5+phase4],[0,ymax],'--',color='0.7')
    ax3.text(0.5+phase4+off, frac*ymax, "$\phi'_4$",fontsize=16)

    plt.savefig('eclipse-geometry3.png',bbox_inches='tight')


