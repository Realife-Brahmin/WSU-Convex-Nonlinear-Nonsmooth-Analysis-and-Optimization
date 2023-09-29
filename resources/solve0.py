#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script for solving the Rosenbrock Test Function
Version: Aug  10, 2023
Author: Tom Asaki
"""

import numpy as np
import optimize as opt

from objective import rosenbrock as obj

#x=np.array([[.1234],[.2345],[-.333],[-1.123],[0],[0],[0],[0],[0]])
x=np.random.randn(12,1)
p=[10]
        
alg=dict(method     = 'GradientDescent',
         maxiter    = 2000,
         ngtol      = 0,
         dftol      = 0,
         dxtol      = 0,
         Lambda     = 1,
         Lambdamax  = 100,
         linesearch = 'Armijo',
         c1         = 0.0001,
         c2         = 0.9,
         progress   = 100
         )   

res=opt.minimize(obj,x,p,alg)

opt.ShowResults(res)
