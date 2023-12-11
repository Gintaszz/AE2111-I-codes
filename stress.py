from zenvmer import shear_pls, moment_pls, torque_pls, bounds
from main import I, c
from data import design, C_R, topweb, bottomweb, centroid, airfoilfunc_top, airfoilfunc_bottom
import numpy as np


def stress(y, designindex, skin=False):
    #y= np.linspace(bounds[0], bounds[1], 150)
    chord = c(y)

    #bending stress
    if skin:
        max_z = max([abs(airfoilfunc_top(x)*chord-centroid(chord, designindex)) for x in np.linspace(0, 1, 100)] +
                    [abs(airfoilfunc_bottom(x)*chord-centroid(chord, designindex)) for x in np.linspace(0, 1, 100)])
    else:
        max_z = max(abs(topweb(0, chord, designindex)-centroid(chord, designindex)), abs(bottomweb(0, chord, designindex)-centroid(chord, designindex)))

    normal_stress = np.abs(moment_pls(y) * max_z / I(y, designindex)) #Absoluite value for now !!!!!!!!!!!!!!!! must be changed
    #torque shear stress
    # \tau = q/t = T/(2tA_m)
    #Area of trapezium
    totalarea = (design['front spar to root chord'][designindex]/C_R + design['back spar to root chord'][designindex]/C_R)*chord/2 * design['spar distance x'][designindex]*chord 
    shearflow = torque_pls(y)/(2*totalarea)
    shear_stress_web = shearflow/design['t web'][designindex]
    #Max shear always going to be at minimum thickness --> web (top and bottom)
    shear_stress = shear_stress_web

    #mohrs circle magic, convert into equivalent normal stress
    mohrs_radius = np.sqrt((normal_stress/2)**2 + shear_stress**2)
    equ_normal = mohrs_radius + normal_stress/2

    return equ_normal