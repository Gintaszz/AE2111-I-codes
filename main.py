from data import design, moment, torque, SPAN, SEMISPAN, C_R, TAPER, DIHEDRAL, E, G, \
                    c, centroid, topweb, bottomweb, airfoilfunc_top, airfoilfunc_bottom

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

def n_stringers(x, designindex): 
    #Interpolation of number of stringers list and span.
    return sp.interpolate.interp1d(design['span list'][designindex], design['list stringers'][designindex], kind="previous",fill_value="extrapolate")(x)


def M(x): #From WP4.1
    return moment(x)


def T(x): #From WP4.1
    return torque(x)


def I(y, designindex): #Moment of inertia calculation function (Ixx)
    chord = c(y)
    spardist = design['spar distance x'][designindex] * chord
    frontsparlength = chord * design['front spar to root chord'][designindex]/C_R
    backsparlength = chord * design['back spar to root chord'][designindex]/C_R
    t_spar = design['t spar'][designindex]
    t_web = design['t web'][designindex]
    centroidlocal = centroid(chord, designindex)

    centroidfrontspar = (airfoilfunc_top(design['front spar x'][designindex]) + 
                           airfoilfunc_bottom(design['front spar x'][designindex])) * chord/2
    
    centroidbackspar = (airfoilfunc_top(design['back spar x'][designindex]) + 
                          airfoilfunc_bottom(design['back spar x'][designindex])) * chord/2
    
    topwebdistsquare = sp.integrate.quad(lambda x: (topweb(x, chord, designindex)-centroidlocal)**2, 0, spardist)[0]/(spardist)
    bottomwebdistsquare = sp.integrate.quad(lambda x: (bottomweb(x, chord, designindex)-centroidlocal)**2, 0, spardist)[0]/(spardist)
    # ASSUMPTION stringers evenly placed
    # I_stringers = area * number * mean of distaces squared (for both top and bottom web)
    I_stringers = design['area stringer'][designindex] * n_stringers(y, designindex)/2 * (topwebdistsquare + bottomwebdistsquare)

    I_frontspar = 1/12 * (frontsparlength)**3 * t_spar + \
                (frontsparlength) * t_spar * (centroidfrontspar - centroidlocal)**2
    I_backspar = 1/12 * (backsparlength)**3 * t_spar + \
                (backsparlength) * t_spar * (centroidbackspar - centroidlocal)**2

    betta_top = np.arctan(topweb(1, 0, designindex))
    betta_bottom = np.arctan(bottomweb(1, 0, designindex))
    
    length_top = spardist / np.cos(betta_top)
    length_bottom = spardist / np.cos(betta_bottom)

    I_topweb = length_top * t_web * ((topweb(0, chord, designindex)+topweb(spardist, chord, designindex))/2 - centroidlocal)**2
    I_bottomweb = length_bottom * t_web * ((bottomweb(0, chord, designindex)+bottomweb(spardist, chord, designindex))/2 - centroidlocal)**2

    return I_stringers + I_frontspar + I_backspar + I_topweb + I_bottomweb

#-------------------------------------------------------------------input of Gintas--------------------------------------------------------------------------
def Izz(y, designindex):


    chord = c(y)
    spardist = design['spar distance x'][designindex] * chord
    frontsparlength = chord * design['front spar to root chord'][designindex]/C_R
    backsparlength = chord * design['back spar to root chord'][designindex]/C_R
    t_spar = design['t spar'][designindex]
    t_web = design['t web'][designindex]
    cx, cz = centroid(chord, designindex)

    centroidfrontspar = (airfoilfunc_top(design['front spar x'][designindex]) + 
                           airfoilfunc_bottom(design['front spar x'][designindex])) * chord/2
    
    centroidbackspar = (airfoilfunc_top(design['back spar x'][designindex]) + 
                          airfoilfunc_bottom(design['back spar x'][designindex])) * chord/2
    
    #topwebdistsquare = sp.integrate.quad(lambda x: (topweb(x, chord, designindex)-centroidlocal)**2, 0, spardist)[0]/(spardist)
    topwebdistsquare = ((design['back spar x'][designindex]-(spardist/2))-cx)**2
    #bottomwebdistsquare = sp.integrate.quad(lambda x: (bottomweb(x, chord, designindex)-centroidlocal)**2, 0, spardist)[0]/(spardist)
    bottomwebdistsquare = ((design['back spar x'][designindex]-(spardist/2))-cx)**2
    

    # ASSUMPTION stringers evenly placed
    # I_stringers = area * number * mean of distaces squared (for both top and bottom web)
    #I_stringers = design['area stringer'][designindex] * n_stringers(y, designindex)/2 * (topwebdistsquare + bottomwebdistsquare)
    sec_dist = spardist/(n_stringers(y, designindex)/2)
    I_stringers = sum([design['area stringer'][designindex]*(i*sec_dist-((design['back spar x'][designindex]-spardist)-cx))**2 for i in range(n_stringers(y, designindex)/2) if True])

    I_frontspar = 1/12 * (frontsparlength) * t_spar**3 + (frontsparlength) * t_spar * (design['back spar x'][designindex]-spardist - cx)**2
    I_backspar = 1/12 * (backsparlength) * t_spar**3 + (backsparlength) * t_spar * (design['back spar x'][designindex] - cx)**2

    betta_top = np.arctan(topweb(1, 0, designindex))  #a little iffy still
    betta_bottom = np.arctan(bottomweb(1, 0, designindex))
    
    length_top = spardist / np.cos(betta_top)
    length_bottom = spardist / np.cos(betta_bottom)

    #I_topweb = length_top * t_web * ((topweb(0, chord, designindex)+topweb(spardist, chord, designindex))/2 - centroidlocal)**2
    I_topweb = ((length_top**3) * (t_web) *(np.cos(betta_top)**2)/12) + length_top*t_web*((design['back spar x'][designindex]-(spardist/2))-cx)**2
    I_bottomweb = ((length_bottom**3) * (t_web) *(np.cos(betta_bottom)**2)/12) + length_bottom*t_web*((design['back spar x'][designindex]-(spardist/2))-cx)**2
    #I_bottomweb = length_bottom * t_web * ((bottomweb(0, chord, designindex)+bottomweb(spardist, chord, designindex))/2 - centroidlocal)**2

    return I_stringers + I_frontspar + I_backspar + I_topweb + I_bottomweb

def Ixz(y, designindex):

    #Assume the top and bottom webs are symetryc
    return 0

#---------------------------------------------------------------------------------------------------------------------------------------------


def J(x, designindex): #Torsional constant calculation = (4tA^2)/s
    t_spar = design['t spar'][designindex]
    t_web = design['t web'][designindex]
    chord = c(x)

    front_spar_length = design['front spar to root chord'][designindex] * chord / C_R
    back_spar_length = design['back spar to root chord'][designindex] * chord / C_R
    web_length = design['spar distance x'][designindex] * chord
    areatotal = (front_spar_length + back_spar_length) * web_length / 2

    torsionalconst = (4 * areatotal**2)/(front_spar_length/t_spar + back_spar_length/t_spar + 2*web_length/t_web)
    return torsionalconst


if __name__ == "__main__":
    DESIGNNUM = 1 #0, 1, or 2
    #-----------------------Integration-------------------------
    stepsint = 100

    #Calculate deflection using .antiderivative(2) attribute
    integrand_moment = [(-M(x)/(E * I(x, designindex=DESIGNNUM))) for x in np.linspace(0, SEMISPAN, stepsint)]
    v = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_moment).antiderivative(2)

    #Calculate twist using .antiderivative(1) attribute
    integrand_torsion = [(T(x)/(G * J(x, designindex=DESIGNNUM))) for x in np.linspace(0, SEMISPAN, stepsint)]
    phi = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_torsion).antiderivative(1)


    #Easy over the limit determination
    deflection = v(SEMISPAN)
    """
    print(f"Maximum deflection of wing: {deflection}\r\nMaximum allowed deflection: {0.15*SPAN}\r\nDifference: {-deflection + 0.15*SPAN}, {(deflection)/(0.15*SPAN)*100}%")
    if -deflection + 0.15*SPAN < 0:
        print("\033[1;31;40mToo much deflection \033[0m")
    """
    twist = phi(SEMISPAN)
    """
    print(f"\r\nMaximum torsion of wing: {twist}\r\nMaximum allowed angle: {10 * np.pi/180}\r\nDifference: {-twist + 10 * np.pi/180}, {(twist)/(10 * np.pi/180)*100}%")
    if -twist + 10 * np.pi/180 < 0:
        print("\033[1;31;40mToo much twist \033[0m")
    """

    print(I(0, designindex=0), J(0, designindex=0))
    print(I(0, designindex=1), J(0, designindex=1))
    print(I(0, designindex=2), J(0, designindex=2))

    plottingres = 0.01

    #-----------------Plotting Wing Deflection-----------------
    '''
    plt.plot(np.arange(0, SEMISPAN, plottingres), [(i*np.sin(DIHEDRAL)) for i in np.arange(0, SEMISPAN, plottingres)], color = 'black', label = "Undeflected wing")

    integrand_moment = [(-M(x)/(E * I(x, designindex=0))) for x in np.linspace(0, SEMISPAN, stepsint)]
    v = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_moment).antiderivative(2)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [(v(i) + i*np.sin(DIHEDRAL)) for i in np.arange(0, SEMISPAN, plottingres)], color = 'b', label = "Design 1")
    print(v(SEMISPAN))

    integrand_moment = [(-M(x)/(E * I(x, designindex=1))) for x in np.linspace(0, SEMISPAN, stepsint)]
    v = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_moment).antiderivative(2)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [(v(i) + i*np.sin(DIHEDRAL)) for i in np.arange(0, SEMISPAN, plottingres)], color = 'orange', label = "Design 2")
    print(v(SEMISPAN))

    integrand_moment = [(-M(x)/(E * I(x, designindex=2))) for x in np.linspace(0, SEMISPAN, stepsint)]
    v = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_moment).antiderivative(2)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [(v(i) + i*np.sin(DIHEDRAL)) for i in np.arange(0, SEMISPAN, plottingres)], color = 'g', label = "Design 3")
    print(v(SEMISPAN))

    plt.axhline(y = 0.15 * SPAN + SEMISPAN * np.sin(DIHEDRAL), color = 'r', linestyle = 'dashed', label = "Limit")
    plt.title("Deflected Wing Front View")
    plt.ylabel("Height [m]")
    plt.xlabel("Span [m]")
    plt.legend()
    plt.show()
    '''

    #-----------------Plotting Wing Twist-----------------
    '''
    integrand_torsion = [(T(x)/(G * J(x, designindex=0))) for x in np.linspace(0, SEMISPAN, stepsint)]
    phi = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_torsion).antiderivative(1)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [phi(i) for i in np.arange(0, SEMISPAN, plottingres)], label = "Design 1")
    print(phi(SEMISPAN))

    integrand_torsion = [(T(x)/(G * J(x, designindex=1))) for x in np.linspace(0, SEMISPAN, stepsint)]
    phi = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_torsion).antiderivative(1)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [phi(i) for i in np.arange(0, SEMISPAN, plottingres)], label = "Design 2")
    print(phi(SEMISPAN))

    integrand_torsion = [(T(x)/(G * J(x, designindex=2))) for x in np.linspace(0, SEMISPAN, stepsint)]
    phi = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, stepsint), integrand_torsion).antiderivative(1)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [phi(i) for i in np.arange(0, SEMISPAN, plottingres)], label = "Design 3")
    print(phi(SEMISPAN))

    plt.axhline(y = -10 * np.pi/180, color = 'r', linestyle = 'dashed', label = "Limit")
    plt.ylabel("Torsion angle [rad]")
    plt.xlabel("Span [m]")
    plt.legend()
    plt.show()
    '''

    #Plotting MOI
    '''
    plt.plot(np.arange(0, SEMISPAN, plottingres), [I(i, 0) for i in np.arange(0, SEMISPAN, plottingres)], label = f'Design 1')
    plt.plot(np.arange(0, SEMISPAN, plottingres), [I(i, 1) for i in np.arange(0, SEMISPAN, plottingres)], label = f'Design 2')
    plt.plot(np.arange(0, SEMISPAN, plottingres), [I(i, 2) for i in np.arange(0, SEMISPAN, plottingres)], label = f'Design 3')
    plt.ylabel("MoI")
    plt.xlabel("Span [m]")
    plt.legend()

    plt.show()
    '''

    #------------Plot 4 graphs for chosen design------------------

    plt.subplot(221)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [(i*np.sin(DIHEDRAL)) for i in np.arange(0, SEMISPAN, plottingres)], color = 'g')
    plt.plot(np.arange(0, SEMISPAN, plottingres), [(v(i) + i*np.sin(DIHEDRAL)) for i in np.arange(0, SEMISPAN, plottingres)])
    plt.axhline(y = 0.15 * SPAN + SEMISPAN * np.sin(DIHEDRAL), color = 'r', linestyle = 'dashed')
    plt.ylabel("Front Profile of wing")
    plt.ylim(-0.2, max(0.15 * SPAN + SEMISPAN * np.sin(DIHEDRAL), deflection + SEMISPAN*np.sin(DIHEDRAL)) + 0.2)
    plt.xlabel("Span [m]")

    plt.subplot(222)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [phi(i) for i in np.arange(0, SEMISPAN, plottingres)])
    plt.axhline(y = -10 * np.pi/180, color = 'r', linestyle = 'dashed')
    plt.ylabel("Torsion angle [rad]")
    plt.xlabel("Span [m]")

    plt.subplot(223)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [I(i, designindex=DESIGNNUM) for i in np.arange(0, SEMISPAN, plottingres)])
    plt.ylabel("MOI")
    plt.xlabel("Span [m]")

    plt.subplot(224)
    plt.plot(np.arange(0, SEMISPAN, plottingres), [J(i, designindex=DESIGNNUM) for i in np.arange(0, SEMISPAN, plottingres)])
    plt.ylabel("Torsional constant")
    plt.xlabel("Span [m]")

    plt.subplots_adjust(hspace=0.3)
    plt.show()