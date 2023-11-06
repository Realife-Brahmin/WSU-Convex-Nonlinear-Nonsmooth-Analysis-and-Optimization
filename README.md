# MATH-564-Convex-and-Nonlinear-Optimization
Convex and Nonlinear Optimization | Fall 2023 | Washington State University | Prof. Thomas Asaki

This repo contains a solver in julia , and different types of optimization problems to solve for.

## Projects as part of the course:
- Parameter Estimation (Estimating DampedSHM parameters)
- Functional Estimation (Estimating a function which minimizes drag force)
- Signal Denosing (Fitting a curve to data, choosing between abiding smoothness of curve or accuracy of fitting)
- Neural Network Training (Detecting Liver Disease in a dataset of patients based on 10 features)
- ???

## LineSearch/TrustRegion Methods Used:
- Conjugate Gradient Descent
- Gradient Descent
- Quasi Newton BFGS
- SR1 based Trust Region Method

## Step-selection Algorithms Used:
- Strong Wolfe + Bisection Interpolation
- Armijo + Backtracking (support removed, always use Strong Wolfe henceforth)
