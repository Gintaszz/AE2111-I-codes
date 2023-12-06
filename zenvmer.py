import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate
from config import *

from mpl_toolkits.axisartist.axislines import AxesZero



#ze dynamic optimizition for reduzed computatioon
zed = {'V': {}, 'M': {}, 'T': {}}
#print(wght)




def chord(x):
  b = wingspan
  lambd = tapper
  c = Cr * (1 - lambd) * 2 * x
  c /= b
  c -= Cr
  return -c


def sweep(c):  #c-chordwise location of sweep(0-1)
  L025 = np.deg2rad(18.6)
  if c > 0.249 and c < 0.251:
    return L025
  ar = 13.5
  lambd = 0.335
  return float(
      np.arctan(
          np.tan(L025) - (4 / ar) * (((c - 0.25) * (1 - lambd)) /
                                     (1 + lambd))))

def coffs():#a0y^2 + a1y + a2
  if('coffs' in zed):
    return zed['coffs'][0],zed['coffs'][1],zed['coffs'][2]
  wrho = (W_wing[0]+W_fuelsys[0]+maxFuel*wght[0])/(wingVolume) #wing density
  beta2=(wrho*(pos_airfoil_thickness+3)/6)*(airfoil_thickness)*(np.cos(sweep(0.25))**2)*(chord(0)**2)
  zed['coffs']=(4*beta2*((1-tapper)**2)/(wingspan**2),-4*beta2*((1-tapper))/(wingspan),beta2)
  return zed['coffs'][0],zed['coffs'][1],zed['coffs'][2]

def wweight(x,fusegroup=True):  #[N] excludes mlg,
  #W_fusegroup = W_vtail[0] + W_htail[0] + W_fuselage[0] + W_nlg[0] + W_engines[0] + W_hydraulics[0] + W_flight_controls[0] + W_electrical[0] + W_avionics[0] + W_ac_and_ai[0]  + W_furnishing[0]+wght[1]*maxPayload

  #make shit symetric around x=0
  x = np.abs(x)
  a0,a1,a2 = coffs()
  if type(x) == type(1.0):

    #contribution of the wing
    monstrocity = (a0 * x**2) + (a1 *x) + a2  #ze wing
    if(fusegroup):
      #contribution of the fuselage
      if ((x > -1.4) and (x < 1.4)):
        monstrocity += (W_fusegroup / diameter)
        return monstrocity * g
  else:
    monstrocity = np.zeros(x.shape)
    #contribution of the wing
    monstrocity[(x < 14.075) * (x > -14.075)] = (a0 * x**2) + (a1 *x) + a2  #ze wing
    if fusegroup:
      #contribution of the fuselage
      monstrocity[(x > -1.4) * (x < 1.4)] += (W_fusegroup / diameter)
  return monstrocity * g


#----------------------------------------------------------------------------------
#READING DATA

a0 = np.array(
    np.genfromtxt('Data.txt',
                  skip_header=21,
                  skip_footer=1029,
                  usecols=(0, 3, 5, 7)))  #[[y1,l1,d1,m1],.. .]
a10=np.array(
    np.genfromtxt('Data10.txt',
                  skip_header=21,
                  skip_footer=1029,
                  usecols=(0, 3, 5, 7)))  #[[y1,l1,d1,m1],.. .]

bounds = (0, a0[-1, 0])


#Maybe make readable
C_L0=0.280980
C_L10=1.214283
C_D0=0.001854
C_D10=0.034000
C_M0=-0.400971
C_M10=-1.54054

alpha = np.arcsin(((C_Ld-C_L0)/(C_L10-C_L0))*np.sin(np.deg2rad(10)))

C_Ldy=a0[:, 1]+((C_Ld-C_L0)/(C_L10-C_L0))*(a10[:, 1]-a0[:,1])
#C_Ddy=a0[:, 2]+((C_Dd-C_D0)/(C_D10-C_D0))*(a10[:, 2]-a0[:,2])
C_Ddy=a0[:, 2]+((C_Ld-C_L0)/(C_L10-C_L0))*(a10[:, 2]-a0[:,2])
#C_Mdy=a0[:, 3]+((C_Md-C_M0)/(C_M10-C_M0))*(a10[:, 3]-a0[:,3])
#C_Md=((C_M10-C_M0)/(np.deg2rad(10)))*alpha+C_M0
#C_Mdy=a0[:, 3]+((C_Md-C_M0)/(C_M10-C_M0))*(a10[:, 3]-a0[:,3])
C_Mdy=a0[:, 3]+((C_Ld-C_L0)/(C_L10-C_L0))*(a10[:, 3]-a0[:,3])






#making coefficients into forces (very wroooooong)
C_Ldy *= q * chord(np.abs(a0[:,0]))
C_Ddy *= q * chord(np.abs(a0[:,0]))
C_Mdy *= q * (chord(np.abs(a0[:, 0]))**2)

#INTERPOLATION
l = sp.interpolate.interp1d(list(a0[:, 0]),
                            list(C_Ldy),
                            kind='cubic',
                            fill_value="extrapolate")
dr = sp.interpolate.interp1d(a0[:, 0],
                            C_Ddy,
                            kind='cubic',
                            fill_value="extrapolate")
m = sp.interpolate.interp1d(a0[:, 0],
                            C_Mdy,
                            kind='cubic',
                            fill_value="extrapolate")

alpha = np.arcsin(((C_Ld-C_L0)*np.sin(np.deg2rad(10)))/(C_L10-C_L0))
#GRAPHING GENERAL
#fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)
fig = plt.figure(figsize=(6.4,5.2))
ax1,ax2,ax3 = fig.add_subplot(311,axes_class=AxesZero),fig.add_subplot(312,axes_class=AxesZero),fig.add_subplot(313,axes_class=AxesZero)

for ax in [ax1,ax2,ax3]:
  for direction in ["xzero", "yzero"]:
    # adds arrows at the ends of each axis
    ax.axis[direction].set_axisline_style("-|>")
    #print(ax)

    # adds X and Y-axis from the origin
    ax.axis[direction].set_visible(True)

  for direction in ["left", "right", "bottom", "top"]:
    # hides borders
    ax.axis[direction].set_visible(False)
  ax.grid()


dihedral = np.deg2rad(1.3)
x = np.linspace(bounds[0], bounds[1], 150)



def axial(x):
  return (dr(x) * np.cos(alpha) - l(x) * np.sin(alpha))#+thrust


def normal(x):
  return l(x) * np.cos(alpha) + dr(x) * np.sin(alpha)



#-------------------------------------------------------------------------------------------------
#Shear Normal (Philosophical considerations)
def shearN(x):
  if not x in zed['V']:
    estimateN, errorg = sp.integrate.quad(lambda x: -(normal(x) - wweight(x) * np.cos(alpha)) * np.cos(dihedral),x, bounds[1])
    if x<=2.5:
        estimateN-=W_mlg[0]*g/2
    zed['V'][x] = estimateN
  return zed['V'][x]

#drag, lift --mathemagic--> normal force + weight(normal component) --cos(dihedral)--> force normal to the local surface
ax1.plot(x, [shearN(i)*n/1000 for i in x if True], 'b-')
v = sp.interpolate.interp1d(list(zed['V'].keys()),
                            list(zed['V'].values()),
                            kind='cubic',
                            fill_value="extrapolate")
#ax1.title.set_text('Normal shear stress')
ax1.set_xlabel(r'y, [$m$]')
ax1.set_ylabel(r'Internal shear, [$kN$]')
ax1.set_xlim(0,14.2)

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#Moments
def momentN(x):  
  #derivative of moment = magnitude of the shear forces
  if not x in zed['M']:
    estimateM, errorg = sp.integrate.quad(lambda x: v(x), x, bounds[1])
    zed['M'][x] = estimateM
  return zed['M'][x]


ax2.plot(x, [momentN(i)*n/1000 for i in x if True], 'b-')
#ax2.title.set_text('Moment diagram')
ax2.set_xlabel(r'y, [$m$]')
ax2.set_ylabel(r'Internal moment, [$kNm$]')
ax2.set_xlim(0,14.2)







#-------------------------------------------------------------------------------------------------
#Torque diagram:
#-------------------------------------------------------------------------------------------------

#estimateD, errorg = sp.integrate.quad(lambda x: axial(x),0, bounds[1])
#print(estimateD)
#thrust = -estimateD
thrust = -31.32*1000
def d_c(x):  # distance of the shear center to some point (x,y) (x = 0 at root quarter chord)
  posx=0.445687 #Goood???
  return (posx-0.25)*chord(x)
def d(x,y):# distance of the shear center to some point (x,y) (x = 0 at root quarter chord)
  posx=0.445687 #Goood???
  a=-(np.tan(-sweep(1-posx)))  
  return (a*x+y)/(((a**2)+1)**0.5)
def dz(): #for distance between engine and shear center
  posz=0.01915
  engineDistance=1.6#m
  return engineDistance
def distributedForce(dist,y):
  return (dist[0] if (dist[1][1]>=y and dist[2][1]<=y) else 0)*d(dist[2][1],y) #wrong but works (make the d(of point in between the two points, now its assumed, two points are on the same x pos))

def torque(x):  #(no point loads or point moments)
  #So things are either distributed or points
  #---------------------------
  #here we group distributions:
  distrits = []
  #distrits+=(magnitude,(start x, start y),(stop x, stop y))

  #fuselagegroup weight distributed along at cg between -r and r
  distrits+=[(-np.cos(alpha) * np.cos(dihedral)*g*W_fusegroup/diameter,(W_fuselage[1],(-diameter/2)),(W_fuselage[1],diameter/2))]
  #---------------------------
  #here we group point weights
  pnts=[]
  #pnts+=(magnitude,(x,y))
  #pnts+=[(-(W_mlg[0]*g/2),(W_mlg[1],W_mlg[2]))]
  pnts+=[(-np.cos(alpha) * np.cos(dihedral)*g*W_mlg[0]/2,(W_mlg[1],W_mlg[2]))]
  
  #---------------------------
  if not x in zed['T']:
    #the integral:
    estimateT, errorg = sp.integrate.quad(lambda x: (normal(x) - wweight(x,False) * np.cos(alpha)) * np.cos(dihedral) * d_c(x) + sum([distributedForce(i,x) for i in distrits if True]) + m(x),
                                          x, bounds[1])
    #the point weights:
    for pnt in pnts:
      estimateT+= pnt[0]*d(pnt[1][0],pnt[1][1]) if(x<=pnt[1][1]) else 0
    estimateT+= thrust*dz() if(x<=(diameter/2)) else 0
    zed['T'][x] = estimateT
  return zed['T'][x]

def moment_pls():#returns interpolated moment
  return sp.interpolate.interp1d(x,[momentN(i)*n for i in x if True], kind='cubic', fill_value="extrapolate")
  
def torque_pls():#returns interpolated torque
  return sp.interpolate.interp1d(x,[torque(i)*n for i in x if True], kind='cubic', fill_value="extrapolate")

ax3.plot(x, [torque(i)*n/1000 for i in x if True], 'b-')
#ax3.title.set_text('Torque diagram')
ax3.set_ylabel(r"Internal torque, [$kNm$]")
ax3.set_xlabel(r'y, [$m$]')
ax3.set_xlim(0,14.2)

'''
[torque(i)*n/1000 for i in x if True]
ax3.plot(x, [-wweight(i)*np.cos(alpha)*n/1000 for i in x if True],)
ax3.plot(x, [normal(i)*n/1000 for i in x if True],)
ax3.arrow(2.5,-wweight(2.5)*np.cos(alpha)*n/1000,0,-(W_mlg[0]/2)*g/1000)
ax3.set_ylabel(r"Loading, [$kN$]")
ax3.set_xlabel(r'y, [$m$]')
ax3.set_xlim(0,14.2)


est, err = sp.integrate.quad(lambda x: wweight(x), 0, bounds[1])

print(f' for weight we have {(est*2)+(W_mlg[0]*g)}')
est, err = sp.integrate.quad(lambda x: l(x), 0, bounds[1])
print(f' for lift we have {(est)*2}')
print(f' or maybe we have L={C_Ld*q*S}, W={weight(wght)*g}')

f = open("demofile2.txt", "a")
f.write(f'{zed['M'].values()}')
f.write('\n')
f.write(f'{zed['T'].values()}')
f.close()

print(f"maximum shear is {min(zed['V'].values())*n}")
print(f"maximum moment is {min(zed['M'].values())*n}")
print(f"maximum torque is {min(zed['T'].values())*n}")
'''
#plt.subplots_adjust(hspace=0.5)

plt.subplots_adjust(bottom=0.07,top=0.97,hspace=0.5)

if __name__ == '__main__':
  plt.show()

#no graphing just return interpolated shit
