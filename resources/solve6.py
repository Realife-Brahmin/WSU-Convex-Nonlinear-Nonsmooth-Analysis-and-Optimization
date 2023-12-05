#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script for solving the Minimium Time Path Problem
Version: November 9, 2023
Author: Tom Asaki
"""

import numpy as np
import optimize as opt
import pandas as pd
import matplotlib.pyplot as plt

from objective import pathtime as obj

N=12  # This is the order of the fit (2N decision variables)

# read in the velocity data array defined on 
# [0,1]x[0,1] and set the path end points
v=pd.read_csv('SpeedData.csv',header=None).to_numpy()
my,mx=v.shape
A=(.05,.05)
B=(.95,.95)

alg=dict(obj        = obj,
         x0         = 0.2*np.random.randn(2*N,1),
         params     = (v,A,B),
         method     = 'BFGS',
         maxiter    = 999,
         progress   = 10,
         ngtol      = 1E-8,
         dftol      = 1E-8,
         dxtol      = 1E-8,
         Lambda     = 1.,
         Lambdamax  = 100.,
         linesearch = 'StrongWolfe',
         c1         = 0.001,
         c2         = 0.9,
         m          = 10,
         maxcond    = 1000,
         )   

res=opt.minimize(alg)

######################################################################
# plot the optimal path superimposed on the velocity image
smp         = 1000      # number of points defining the path
FigDPI      = 256       # figure dpi (effects scale)
FigSize     = (8,6)     # fiugure size   
ColorMap    = 'jet'     # velocity colormap
LineColor   = 'white'   # path plot color
LineWidth   = 1         # path line width
PointSize   = 16        # size of path endpoints

r=np.linspace(0,1,smp)
xx=(1-r)*A[0]+r*B[0]
yy=(1-r)*A[1]+r*B[1]
for k in range(N):
    s=np.sin((k+1)*np.pi*r)
    xx+=res['x'][k,-1]*s
    yy+=res['x'][k+N,-1]*s
xxr=xx*(mx-1)
yyr=yy*(my-1)

fig = plt.figure(dpi=FigDPI,figsize=FigSize)
ax=fig.add_subplot()
vim=ax.imshow(v,cmap=ColorMap)
plt.colorbar(vim,orientation='vertical')
ax.plot(yyr,xxr,color=LineColor,linewidth=LineWidth)
ax.scatter(A[1]*mx,A[0]*my,PointSize,LineColor)
ax.scatter(B[1]*mx,B[0]*my,PointSize,LineColor)
plt.xticks([])
plt.yticks([])
plt.show()



