# Maximum allowable stress as a function of span
# Web, Skin, Column, Global Buckling
import numpy as np
from data import design, C_R, topweb, bottomweb, centroid, c, SEMISPAN, E
import matplotlib.pyplot as plt

def skinmaxstress(y, designindex):
    chord = c(y)
    frontsparlength = design['front spar to root chord'][designindex]/C_R * chord
    backsparlength = design['back spar to root chord'][designindex]/C_R * chord
    spardist = design['spar distance x'][designindex] * chord

    k_c = 7
    t = design['t web'][designindex]
    b = (frontsparlength + backsparlength)/2
    poission = 1/3
    sigma = np.pi**2 * k_c * E * (t/b)**2 / (12(1-poission**2))

    return sigma

plotlist = [skinmaxstress(y, 0) for y in np.arange(0, SEMISPAN, 0.1)]
plt.plot(np.arange(0, SEMISPAN, 0.1), plotlist)
plt.show()