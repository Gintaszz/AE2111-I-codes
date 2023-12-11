from stress import stress
from config import *
from mpl_toolkits.axisartist.axislines import AxesZero
import matplotlib.pyplot as plt
from zenvmer import bounds

sigmaY = 276000000.0
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
  #ax.set_yscale('log')
#x = np.linspace(bounds[0], bounds[1], 150)
x = np.linspace(bounds[0], 7, 150)

ax1.plot(x,[sigmaY/stress(i,0) for i in x if True])
ax2.plot(x,[sigmaY/stress(i,1) for i in x if True])
ax3.plot(x,[sigmaY/stress(i,2) for i in x if True])

if __name__ == '__main__':
  plt.subplots_adjust(bottom=0.07,top=0.97,hspace=0.5)
  plt.show()

