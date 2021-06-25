#!/usr/bin/env python

RES = {
    'data01.txt' : (87.5443, 0.09505, 0.21655, 10.0094, 3.5918),
    'data02.txt' : (88.7398, 0.16077, 0.24197, 12.4736, 2.9267),
    'data03.txt' : (81.8334, 0.06411, 0.29096, 5.6060, 2.2509),
    'data04.txt' : (89.4872, 0.12765, 0.13940, 10.1068,  1.7579),
    'data05.txt' : (89.5093, 0.15045, 0.16198, 8.9476, 2.6993),
    'data06.txt' : (87.2902, 0.17918, 0.29453, 5.3872, 2.1097),
    'data07.txt' : (84.1118, 0.12861, 0.29129, 12.4202, 1.7936),
    'data08.txt' : (89.4056, 0.17624, 0.28962, 6.7535, 2.7553),
    'data09.txt' : (85.1774, 0.10259, 0.23282, 14.2960, 2.8617),
    'data10.txt' : (89.7266, 0.19543, 0.23109, 11.9497, 1.8421),
}

MOD = {
    'data01.txt' : (87.6162 , 0.0950, 0.2162, 10.041, 3.601),
    'data02.txt' : (88.8047 , 0.1610, 0.2419, 12.449, 2.923),
    'data03.txt' : (83.0359 , 0.0621, 0.2819, 5.978, 2.397),
    'data04.txt' : (89.4595 , 0.1279, 0.1393, 10.089, 1.744),
    'data05.txt' : (89.5123 , 0.1514, 0.1616, 8.872, 2.675),
    'data06.txt' : (87.2253 , 0.1791, 0.2944, 5.389, 2.113),
    'data07.txt' : (84.1395 , 0.1289, 0.2909, 12.373, 1.797),
    'data08.txt' : (89.6150 , 0.1761, 0.2895, 6.766, 2.756),
    'data09.txt' : (85.2083 , 0.1023, 0.2325, 14.367, 2.870),
    'data10.txt' : (89.8887 , 0.1954, 0.2307, 11.962, 1.838),
}

TITLES = (
    '$i$', '$r_1$', '$r_2$', '$s_1$', '$s_2$'
)

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    fig,axs = plt.subplots(1,5,figsize=(10,2.2))
    
    for n,ax in enumerate(axs):
        xs, ys =[], []
        for key in MOD:
            xs.append(MOD[key][n])
            ys.append(RES[key][n])
        xs = np.array(xs)
        ys = np.array(ys)
        lo = min(xs.min(), ys.min())
        hi = max(xs.max(), ys.max())
        rng = hi-lo
        lo -= 0.1*rng
        hi += 0.1*rng
        ax.plot([lo,hi],[lo,hi],'--',color='0.7')
        ax.plot(xs,ys,'.b')
        ax.set_title(TITLES[n])
        if n == 0:
            ax.set_ylabel('Measured value')
        ax.set_xlabel('Model value')
    plt.tight_layout()
    plt.savefig('compare.png')
            
