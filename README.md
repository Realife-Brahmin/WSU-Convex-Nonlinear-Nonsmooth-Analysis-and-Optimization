# MATH 564: Convex and Nonlinear Optimization and MATH 565: Nonsmooth Analysis and Optimization
## Prof. Thomas Asaki | Fall 2023 and Spring 2024
<img src = "https://s3.wp.wsu.edu/uploads/sites/3144/2023/03/fall-preview-photo.jpg" height = 250pt> <img src = "https://s3.wp.wsu.edu/uploads/sites/227/2022/10/photo_serv_3668_large.jpg" height = 250pt>
 <!--  
 <img src = "https://user-images.githubusercontent.com/24756405/237030072-5b15d383-9fec-4af8-bf5d-572b6db31e37.png" width = 35% height = 30%> <img src = "https://user-images.githubusercontent.com/24756405/237032599-14edd7bc-5b4c-4f0d-84a3-3ff7e9f963ec.png" width = 50% height = 100%>
-->
## Julia implementations for the two courses at Washington State University, Pullman.

This repo contains a solver in julia , and different types of optimization problems to solve for.

## Projects as part of the two courses:
- Parameter Estimation (Estimating DampedSHM parameters)
- Functional Estimation (Estimating a function which minimizes drag force)
- Signal Denosing (Fitting a curve to data, choosing between abiding smoothness of curve or accuracy of fitting)
- Neural Network Training (Detecting Liver Disease in a dataset of patients based on 10 features)
- Minimum Time Trajectory (for a 2D Speed Matrix, and given points A and B on the plane, compute the trajectory minimizing the time of travel)
- Final Project (Computing the optimal placement of a receiver, given the location of several transmitters)
- Derivative Free Optimization Methods (implementation of Genetic Algorithm and Nelder Mead Simplex Method)
- Equality Constrained Quadratic Programming (detecting location of fire based on some boundary constraints which can be expressed as an ECQP)
  
## LineSearch/TrustRegion Methods Used:
- Gradient Descent
- Quasi Newton BFGS
- Projected Gradient Conjugate Gradient Method
- Note: Conjugate Gradient Descent Method hasn't been supported in a while, and has issues.
- Note: SR1 based Trust Region Method implementation has been indefinitely paused, in favour of more urgent tasks.

## Heuristic Methods Available:
- Genetic Algorithm
- Nelder Mead Simplex Algorithm

## Space Sampling Methods Available:
- Halton Sequence Sampling
- Latin Hypercube Sampling
  
## Step-selection Algorithms Used:
- Strong Wolfe + Bisection Interpolation
- Armijo + Backtracking (support removed, always use Strong Wolfe henceforth

<img src = "https://user-images.githubusercontent.com/24756405/237028282-16ffc864-98f8-4c6f-a663-dfb04d191623.png" width = 15% height = 15%>

