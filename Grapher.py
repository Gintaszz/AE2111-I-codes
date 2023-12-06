import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate


#READING DATA
print(np.genfromtxt("Data.txt",skip_header=5,skip_footer=1095)) 
for i in range(10000):
    alpha = i         
#alpha = np.deg2rad(np.genfromtxt("Data.txt",skip_header=5,skip_footer=1094,usecols=(2))) #make alpha also read from file
print(alpha)
a = np.array(
    np.genfromtxt('Data.txt',
                  skip_header=21,
                  skip_footer=1040,
                  usecols=(0, 3, 5, 7)))

print('read succesfully')


#INTERPOLATION
l = sp.interpolate.interp1d(list(a[:, 0]),
                            list(a[:, 1]),
                            kind='cubic',
                            fill_value="extrapolate")
d = sp.interpolate.interp1d(a[:, 0],
                            a[:, 2],
                            kind='cubic',
                            fill_value="extrapolate")
m = sp.interpolate.interp1d(a[:, 0],
                            a[:, 3],
                            kind='cubic',
                            fill_value="extrapolate")




#GRAPHING GENERAL
fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)

#Shear Normal
x=np.linspace(0,a[-1,0],1000)
def shearN(x):
  estimateN,errorg = sp.integrate.quad(lambda x: l(x)*np.cos(alpha),a[0,0],x)
  return estimateN
y=[]
for i in x:
  y.append(shearN(i))
ax1.plot(x,y,'o')
ax1.title.set_text('Normal shear stress')

#Shear Axial
def shearA(x):
  estimateA,errorg = sp.integrate.quad(lambda x: d(x)*np.cos(alpha),a[0,0],x)
  return estimateA
y=[]
for i in x:
  y.append(shearA(i))

ax0.plot(x,y,'o')
ax0.title.set_text('Axial shear stress')



plt.subplots_adjust(hspace=0.5)
plt.show()