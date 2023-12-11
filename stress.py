from zenvmer import shear_pls, moment_pls, torque_pls, bounds
from main import I, c
from data import design, C_R, topweb, bottomweb, centroid, airfoilfunc_top, airfoilfunc_bottom, SEMISPAN
import numpy as np
import matplotlib.pyplot as plt
from functools import cache

@cache
def stress(y, designindex):
    #y= np.linspace(bounds[0], bounds[1], 150)
    chord = c(y)

    #bending stress

    max_z = max([abs(airfoilfunc_top(x)*chord-centroid(chord, designindex)) for x in np.linspace(0, 1, 100)] +
                [abs(airfoilfunc_bottom(x)*chord-centroid(chord, designindex)) for x in np.linspace(0, 1, 100)])

    normal_stress = np.abs(moment_pls(y) * max_z / I(y, designindex)) #Absoluite value for now !!!!!!!!!!!!!!!! must be changed

    
    #torque shear stress
    # \tau = q/t = T/(2tA_m)
    #Area of trapezium
    totalarea = (design['front spar to root chord'][designindex]/C_R + design['back spar to root chord'][designindex]/C_R)*chord/2 * design['spar distance x'][designindex]*chord 
    shearflow_torsion = torque_pls(y)/(2*totalarea)
    shear_stress_web = shearflow_torsion/design['t web'][designindex]
    #Max shear always going to be at minimum thickness --> web (top and bottom)
    #Shear stress due to shear force at corner of wingbox assuming square geometry
    shearflow_shear_corner = -shear_pls(y)/I(y, designindex) * design['t web'][designindex] * design['spar distance x'][designindex] * chord

    shear_stress = shear_stress_web + shearflow_shear_corner

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    shear_stress = 0 #assume no shear in skin
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    #mohrs circle magic, convert into equivalent normal stress
    mohrs_radius = np.sqrt((normal_stress/2)**2 + shear_stress**2)
    equ_normal = mohrs_radius + normal_stress/2

    return equ_normal

if __name__ == "__main__":
    plotlist = [stress(y, 1) for y in np.arange(0, SEMISPAN, 0.1)]
    print(max(plotlist))

    plt.plot(np.arange(0, SEMISPAN, 0.1), plotlist)
    plt.show()