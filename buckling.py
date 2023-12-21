# Maximum allowable stress as a function of span
# Web, Skin, Column, Global Buckling
import numpy as np
from data import design, C_R, c, SEMISPAN, E, tweb
from stress import stress, max_stress_stringer
import matplotlib.pyplot as plt
from functools import cache

@cache
def skinmaxstress(y, designindex, safetyfactor = 1):
    chord = c(y)
    spardist = design['spar distance x'][designindex] * chord

    b = spardist
    
    k_c = 7
    t = tweb(y, designindex) # [m]

    poission = 1/3
    sigma = np.pi**2 * k_c * E * (t/b)**2 / (12*(1-poission**2))

    return sigma/safetyfactor

def webmaxstress(y, designindex, safetyfactor = 1):
    sigma = 0
    return sigma/safetyfactor

def columnmaxstress(y, designindex, safetyfactor = 1):
    t = design['t stringers'][designindex]
    A = design['area stringer'][designindex]
    b = (A + t**2) / (2 * t)

    ybar = 0
    
    I_stringer = 1/12 * b**3 * t + (b/2-ybar)**2*b*t +\
             1/12 * t**3 * (b-t) + (t/2-ybar)**2*(b-t)*t
    
    L_stringer  = design['rib spacing'][designindex] #Rib Spacing
    K = 0.25

    sigma = K * np.pi**2 * E * I_stringer / (L_stringer**2 * A)

    return sigma/safetyfactor

def compressivestress(y, designindex, safetyfactor = 1):
    sigma = 310 * 10**6 # Maximum tensile strength of AL6061-T6
    return sigma/safetyfactor


if __name__ == "__main__":

    for designindex in [0, 1, 2]:
        print(f"Design {designindex}")

        spanarray = np.arange(0, SEMISPAN, 0.1)

        skinmargin = [skinmaxstress(y, designindex)/stress(y, designindex) for y in spanarray]

        webmargin = [webmaxstress(y, designindex)/stress(y, designindex) for y in spanarray]

        sigmamaxcolumn = columnmaxstress(0, designindex)
        columnmargin = [sigmamaxcolumn/max_stress_stringer(y, designindex) for y in spanarray]

        sigmamaxcompressive = compressivestress(0, designindex)
        compressivemargin = [sigmamaxcompressive/stress(y, designindex) for y in spanarray]

        bucklingmargins = {'skin': skinmargin, 'web': webmargin, 'column': columnmargin, 'compressivestress': compressivemargin}

        for mode, modelist in bucklingmargins.items():
            if min(modelist) < 1.5:
                print(f"\033[91mFailed: {mode} buckling \033[0m: {min(modelist)}")
            else:
                print(f"\033[92mSuccess: {mode} buckling \033[0m: {min(modelist)}")
            
            plt.plot(spanarray, modelist, label = f'{mode} buckling')
        
        plt.axhline(1.5, color = 'r', linestyle = ':', label = 'Limit')
        plt.yscale("log")
        plt.ylim(0.5, 10**3)
        plt.ylabel('Margin of Safety [-]', fontsize=16)
        plt.xlabel('Span [m]', fontsize=16)
        plt.legend()
        plt.grid()
    
        plt.show()