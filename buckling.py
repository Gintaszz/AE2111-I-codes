# Maximum allowable stress as a function of span
# Web, Skin, Column, Global Buckling
import numpy as np
from data import design, C_R, topweb, bottomweb, centroid, c, SEMISPAN, E
from main import Izz
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
    t_stringer = 0.001 # 1mm
    A_stringer = design['area stringer'][designindex]
    b_stringer = A_stringer / (2 * t_stringer)
    I_stringer = 4/3 * (t_stringer * b_stringer**3)
    L_stringer = SEMISPAN #Losnges stringer is semispan long
    K = 0.25

    sigma = K * np.pi**2 * E * Izz(y) / (L_stringer* A_stringer**2)
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

    bucklingmargins = {'skin': skinmargin, 'web': webmargin, 'column': columnmargin, 'global': globalmargin}

    for mode, modelist in bucklingmargins.items():
        if min(modelist) < 1:
            print(f"\033[91m Failed: {mode} buckling \033[0m: {min(modelist)}")
        else:
            print(f"\033[92m Success: {mode} buckling \033[0m: {min(modelist)}")
        
        plt.plot(spanarray, modelist)

    plt.ylim(0, 10)
    plt.show()