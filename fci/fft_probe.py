import os
import sys
import numpy as np
import math
#-----------------------------------------------------------
import matplotlib.pyplot as plt
from meshgrid2 import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from invdisttree import Invdisttree
from g_fit3d import GFit
#-----------------------------------------------------------

x = np.linspace(-3.8, 3.8, 200)
y = np.linspace(0.0, 7.34, 200)
z = np.linspace(-3.8, 3.8, 200)

xi,yi,zi = meshgrid2(x,y,z)
XX=np.vstack((xi.flatten(), yi.flatten(), zi.flatten()))

X,F=GFit.read_from_file(0,0,path='/data/users/mklymenko/science/back/programing/results/')
invdisttree = Invdisttree(X.T,F, leafsize=10, stat=1)
AA=invdisttree(XX.T, nnear=130, eps=0, p=1)

a=AA.reshape(xi.shape)
b=np.fft.rfftn(a,[70,70,70])
b=np.fft.fftshift(b)

fxi,fyi = np.meshgrid(np.arange(70),np.arange(70))
fxi,fyi = fxi/(len(x)*(x[3]-x[2])),fyi/(len(y)*(y[3]-y[2]))

fig=plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(fxi,fyi,np.abs(b[:,:,50]))
#ax.plot_surface(xi[:,:,75],yi[:,:,75],a[:,:,75], cmap=cm.jet, linewidth=0.2)
#ax.plot_surface(fxi,fyi,np.abs(b[:,:,0]), cmap=cm.jet)
#plt.contour(fxi, fyi, abs(b[:,:,0]), colors='blue')

plt.imshow(np.abs(b[:,:,0]))
plt.hold(True)
plt.show()
