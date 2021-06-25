#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    phs,fs,fes = np.loadtxt('test_light_curve.txt', unpack=True)

    # measured values
    p3,p4,f0,f1,f2 = 0.015037,0.064282,0.64882,0.29509,0.53987
    
    fig,(axp,axs) = plt.subplots(1,2,figsize=(12,4))

    
    axp.plot([p3,p3],[0.25,0.7],'--',color='0.7')
    axp.plot([p4,p4],[0.25,0.7],'--',color='0.7')
    axp.plot([-0.1,0.1],[f0,f0],'--',color='0.7')
    axp.plot([-0.1,0.1],[f1,f1],'--',color='0.7')

    axp.text(-0.045,f0+0.01,'$f_0$')
    axp.text(-0.045,f1+0.01,'$f_1$')
    axp.text(p3+0.002,0.45,'$\phi_3$')
    axp.text(p4+0.002,0.45,'$\phi_4$')
    axp.errorbar(
        phs,fs,fes,fmt='.b',ecolor='0.7'
    )
    axp.set_xlim(-0.1,0.1)
    axp.set_ylim(0.25,0.7)
    axp.set_xlabel('Orbital phase')
    axp.set_ylabel('Flux')

    axs.plot([0.4,0.6],[f0,f0],'--',color='0.7')
    axs.plot([0.4,0.6],[f2,f2],'--',color='0.7')
    axs.text(0.455,f0+0.01,'$f_0$')
    axs.text(0.455,f2+0.01,'$f_2$')
    axs.errorbar(
        phs,fs,fes,fmt='.b',ecolor='0.7'
    )
    axs.set_xlim(0.4,0.6)
    axs.set_ylim(0.25,0.7)
    axs.set_xlabel('Orbital phase')
    
    plt.tight_layout()
    plt.savefig('annotate.png')
            
