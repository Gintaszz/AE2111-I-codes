# Maximum allowable stress as a function of span
# Web, Skin, Column, Global Buckling
import numpy as np
from data import design, C_R, topweb, bottomweb, centroid, c, SEMISPAN, E
from stress import stress
import matplotlib.pyplot as plt

def skinmaxstress(y, designindex, safetyfactor = 1):
    chord = c(y)
    spardist = design['spar distance x'][designindex] * chord

    b = spardist
    
    k_c = 7
    t = design['t web'][designindex] # [m]

    poission = 1/3
    sigma = np.pi**2 * k_c * E * (t/b)**2 / (12*(1-poission**2))

    return sigma/safetyfactor

def webmaxstress(y, designindex, safetyfactor = 1):
    sigma = 0
    return sigma/safetyfactor

def columnmaxstress(y, designindex, safetyfactor = 1):
    sigma = 0
    return sigma/safetyfactor

def globalmaxstress(y, designindex, safetyfactor = 1):
    sigma = 0
    return sigma/safetyfactor

if __name__ == "__main__":
    designindex = 1
    spanarray = np.arange(0, SEMISPAN, 0.1)

    skinmargin = [skinmaxstress(y, designindex)/stress(y, designindex) for y in spanarray]
    webmargin = [webmaxstress(y, designindex)/stress(y, designindex) for y in spanarray]
    columnmargin = [columnmaxstress(y, designindex)/stress(y, designindex) for y in spanarray]
    globalmargin = [globalmaxstress(y, designindex)/stress(y, designindex) for y in spanarray]

    if min(skinmargin) < 1: print(f"\033[1;31;40m Failed: Skin buckling \033[0m: {min(skinmargin)}")
    if min(webmargin) < 1: print(f"\033[1;31;40m Failed: Web buckling \033[0m: {min(webmargin)}")
    if min(columnmargin) < 1: print(f"\033[1;31;40m Failed: Column buckling \033[0m: {min(columnmargin)}")
    if min(globalmargin) < 1: print(f"\033[1;31;40m Failed: Global buckling \033[0m: {min(globalmargin)}")

    plt.plot(spanarray, skinmargin)
    plt.plot(spanarray, webmargin)
    plt.plot(spanarray, columnmargin)
    plt.plot(spanarray, globalmargin)
    plt.ylim(0, 10)
    plt.show()