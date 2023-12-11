# Maximum allowable stress as a function of span
# Web, Skin, Column, Global Buckling
import numpy as np
from data import design, C_R, topweb, bottomweb, centroid

def skinmaxstress(y):
    chord = c(y)
    
    k_c = (frontsparlength + backsparlength)/2 * spartdist # should be constant
    sigma = np.pi**2 * k_c * E * (t/b)**2 / (12(1-poission**2))
    return