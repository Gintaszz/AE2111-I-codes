# Maximum allowable stress as a function of span
# Web, Skin, Column, Global Buckling
import numpy as np
from data import design, C_R, topweb, bottomweb, centroid, c, SEMISPAN, E
from stress import stress
import matplotlib.pyplot as plt

def skinmaxstress(y, safetyfactor, designindex):
    chord = c(y)
    frontsparlength = design['front spar to root chord'][designindex]/C_R * chord
    backsparlength = design['back spar to root chord'][designindex]/C_R * chord
    spardist = design['spar distance x'][designindex] * chord
    ribspacing = 1

    b = min(ribspacing, spardist)
    
    k_c = 7
    t = 0.015 # [m]

    poission = 1/3
    sigma = np.pi**2 * k_c * E * (t/b)**2 / (12*(1-poission**2))

    return sigma/safetyfactor

plotlist = [skinmaxstress(y, 1, 1)/stress(y, 1, skin=True) for y in np.arange(0, SEMISPAN, 0.1)]
print(min(plotlist))

plt.plot(np.arange(0, SEMISPAN, 0.1), plotlist)
plt.ylim(0, 10)
plt.show()