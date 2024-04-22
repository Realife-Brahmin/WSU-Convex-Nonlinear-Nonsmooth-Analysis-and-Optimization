
"""
        solveForNextPGCGIterate(xk, gk, dk, rk, num, G, Ae, AAT; verbose=false)

Compute the next iterate in a Projected Gradient Conjugate Gradient (PGCG) method.

# Arguments
- `xk::Vector`: Current point in the parameter space.
- `gk::Vector`: Current gradient at point xk.
- `dk::Vector`: Current search direction.
- `rk::Vector`: Current residual, which is the gradient at xk projected onto the feasible set.
- `num::Float64`: The numerator in the calculation of the step size, usually the inner product of rk and gk.
- `G::Matrix`: Hessian matrix of the quadratic objective function or its approximation.
- `Ae::Matrix`: Coefficient matrix for the equality constraints.
- `AAT::Matrix`: Precomputed product of Ae and its transpose (Ae').

# Keyword Arguments
- `verbose::Bool=false`: Flag to enable verbose output, useful for debugging.

# Returns
- `xkp1::Vector`: Next point in the parameter space.
- `gkp1::Vector`: Next gradient.
- `dkp1::Vector`: Next search direction.
- `rkp1::Vector`: Next residual.
- `fevals_1PGCG::Int`: Number of function evaluations made in this step, for bookkeeping.
- `actions_1PGCG::Dict`: A dictionary to hold any additional information on the actions taken, for bookkeeping.

# Description
This function is a single iteration step of the PGCG algorithm, which is a variant of the Conjugate Gradient method adapted for problems with equality constraints. It computes the next iterate by calculating a step size using the formula `alphak = num/den` where `den` is the inner product of the current search direction `dk` and the result of `G*dk`. The algorithm then updates the current point `xk`, the residual `rk`, and the search direction `dk` accordingly. 

# Notes
- It is assumed that `AAT` is non-singular and that it has been computed prior to calling this function to avoid redundant calculations.
- The variables `fevals_1PGCG` and `actions_1PGCG` are returned for further use in bookkeeping and adaptive algorithm adjustments.

# Examples
```julia
# Given a current point xk, gradient gk, search direction dk, and residual rk:
xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG = solveForNextPGCGIterate(
        xk, gk, dk, rk, num, G, Ae, AAT
)
```
"""
function solveForNextPGCGIterate(xk, gk, dk, rk, num, G, Ae, AAT;
        verbose=false)

        den = transpose(dk)*G*dk

        alphak = num/den
        xkp1 = xk + alphak*dk
        rkp1 = rk + alphak*G*dk
        vkp1 = AAT\(Ae*rkp1)
        gkp1 = rkp1 - transpose(Ae)*vkp1
        betak = transpose(rkp1)*gk/num

        dkp1 = -gkp1 + betak*dk

        # adding fevals and actions as oargs for consistency across my opt algos
        fevals_1PGCG = 0
        actions_1PGCG = Dict()

        return xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG
end