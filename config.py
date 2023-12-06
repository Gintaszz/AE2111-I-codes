from ISA import ISA
import numpy as np
#Weights (mass, d from nose, d from root if relevant)

#[kg]
W_wing	=	(1615.867175, 11.2,0)
W_htail = (120.5059257, -21.2+W_wing[1], 0)  #
W_vtail = (310.1311673 , W_wing[1]-20.2,0) #
W_fuselage = (1310.055555,W_wing[1]-8.96)  #revise CG location
W_mlg = (763.8813875, W_wing[1]-13.5,2.5 )#
W_nlg = (77.2604848,W_wing[1]-1.837453962,0)  #
W_engines = (2119.144352,W_wing[1]-18.172,1.65)   #
W_fuelsys	=	(288.5466907,W_wing[1]-11.2,0)
W_flight_controls = (771.8529865,W_fuselage[1])  #
W_hydraulics = (15.78592133,W_fuselage[1])  #
W_electrical = (272.9256765,W_fuselage[1])  #
W_avionics = (604.4872083,W_fuselage[1])  #
W_ac_and_ai = (146.4028914,W_fuselage[1])  #
W_furnishing = (889.2570823,W_fuselage[1])

maxFuel = 5469.4077256#[kg]
maxPayload = 1010#[kg]



def weight(x): #x=wght[kg]
    OEW = W_wing[0]+W_fuelsys[0] +W_vtail[0] + W_htail[0] + W_fuselage[0] + W_nlg[0]+W_mlg[0] + W_engines[0] + W_hydraulics[0] + W_flight_controls[0] + W_electrical[0] + W_avionics[0] + W_ac_and_ai[0] +  W_furnishing[0]
    #print(f"OEW={OEW}")
    return (OEW + x[0]*maxFuel+x[1]*maxPayload)

#Geometry
l_fus = 22.400
diameter = 2.8
MAC=2.8 #revise!!!!!! or scrap
wingVolume=10.0055621 #for some reason this value seems to work, it slightly deviates from the one calculated by Ned
Cr=3.13
pos_airfoil_thickness=0.349
airfoil_thickness=0.15

wingspan = 28.15
tapper = 0.335






#some general shit
S = 58.76
profileDrag = 0.026984543
Aspect = 13.5
oswald = 0.69


g = 9.80665


#------------------------------------------------------------------------------------------------------------
#The load case:
'''
#Dive at sea level **make the diagrams
n=2.64
h = 0#12496.8
V = 262

#Dive at cruise
n=2.64
h = 12496.8
V = 262

#accelerated stall at sea level
n=2.64
h = 0
V = 101

#accelerated stall at sea level
n=2.64
h = 12496.8
V = 208

#inverted cruise
n=-1
h = 12496.8
V = 209

#inverted sea-level
n=-1
h = 0
V = 209
'''
n=1
h = 12496.8
V = 262
wght=(1,1) #(fuel, payload)  (0-1) fraction of maximum that is loaded 
#------------------------------------------------------------------------------------------------------------


p, Temp,rho = ISA(h)
#V= ((1.4 * 287 * Temp)**0.5) * 0.71
q = ((rho * V**2) / 2)

W_fusegroup = W_vtail[0] + W_htail[0] + W_fuselage[0] + W_nlg[0] + W_engines[0] + W_hydraulics[0] + W_flight_controls[0] + W_electrical[0] + W_avionics[0] + W_ac_and_ai[0]  + W_furnishing[0]+wght[1]*maxPayload

C_Ld=weight(wght)*g/(q*S)

#Assume this is right (possibly wrong)
#C_Dd=profileDrag+(C_Ld**2)/(Aspect*np.pi*oswald)

#adding CD0 to CDi from data

#a[:, 2] += profileDrag

#Fucked if I know
#C_Md=-1