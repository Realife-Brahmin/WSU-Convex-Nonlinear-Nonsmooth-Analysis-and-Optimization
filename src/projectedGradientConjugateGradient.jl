"""
        solveForNextPGCGIterate(xk, gk, dk, rk, num, G, A, AAT; verbose=false) -> Tuple

Performs one iteration of the Projected Gradient Conjugate Gradient (PGCG) method for solving ECQP problems, calculating the next values of solution vector, gradient, search direction, and residual.

# Arguments
- `xk`: The current solution vector.
- `gk`: The current gradient of the quadratic objective function with respect to `xk`.
- `dk`: The current search direction.
- `rk`: The current residual vector, which is the gradient of the objective function modified by the projection onto the constraint space.
- `num`: A scalar value used in the calculation of the step size and the update of the search direction.
- `G`: The quadratic coefficient matrix from the objective function.
- `A`: The matrix representing linear equality constraints.
- `AAT`: Precomputed matrix product of `A` and its transpose, used for efficiency in constraint projection.

# Keyword Arguments
- `verbose::Bool=false`: If `true`, enables printing of function-specific messages for debugging or monitoring.

# Returns
- `xkp1`: The next solution vector after the current iteration.
- `gkp1`: The next gradient vector.
- `dkp1`: The next search direction.
- `rkp1`: The next residual vector.
- `fevals_1PGCG`: The number of function evaluations performed in this iteration (currently set to 0 as a placeholder).
- `actions_1PGCG`: A dictionary of actions taken during this iteration (currently an empty dictionary as a placeholder).

This function calculates the step size (`alphak`) for the current search direction (`dk`), updates the solution vector (`xk`), and computes the new residual (`rkp1`) and gradient (`gkp1`). It also computes the next search direction (`dkp1`) based on the updated residual and gradient, incorporating the conjugate gradient method's beta update (`betak`).

# Example
```julia
# Assuming the necessary variables (xk, gk, dk, rk, num, G, A, AAT) are defined
# and correspond to the current state of the ECQP optimization process:

xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG = solveForNextPGCGIterate(
        xk, gk, dk, rk, num, G, A, AAT, verbose=true
)
```
"""
function solveForNextPGCGIterate(xk, gk, dk, rk, num, G, A, AAT;
        verbose=false)

        den = transpose(dk)*G*dk

        alphak = num/den
        xkp1 = xk + alphak*dk
        rkp1 = rk + alphak*G*dk
        vkp1 = AAT\(A*rkp1)
        gkp1 = rkp1 - transpose(A)*vkp1
        betak = transpose(rkp1)*gk/num

        dkp1 = -gkp1 + betak*dk

        # adding fevals and actions as oargs for consistency across my opt algos
        fevals_1PGCG = 0
        actions_1PGCG = Dict()

        return xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG
end