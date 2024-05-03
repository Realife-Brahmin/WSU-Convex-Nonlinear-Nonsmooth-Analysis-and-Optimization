using LinearAlgebra
using Parameters
using Test

include("../helperFunctions.jl")
include("../types.jl")

"""
    checkQualityOfCandidateStep(SR1params, pr, solState, pj, solverState)

Evaluate the quality of a candidate step in a numerical optimization algorithm.

# Arguments
- `SR1params`: Parameters for the SR1 (Symmetric Rank 1) update, typically used in quasi-Newton methods.
  It should contain at least `gkm1` (previous gradient) and `Bkm1` (previous Hessian approximation).
- `pr`: A structure related to the problem being solved, which includes at least the objective function `f`.
- `solState::SolStateType`: Contains the current solution state, including `xk` (current point) and `fk` (function value at `xk`).
- `pj::Vector{Float64}`: The step to be evaluated, proposed from the current point `xk`.
- `solverState::SolverStateType`: Contains state information of the solver, like the count of function evaluations.

# Returns
- A tuple with the quality ratio `ρk` and the updated `solverState`. The ratio `ρk` indicates the quality of the step based on the actual and predicted decrease in the function value.

# Description
This function assesses the quality of a candidate step (`pj`) in an optimization process by comparing the actual decrease in the objective function to the decrease predicted by a quadratic model. It uses the SR1 update parameters and the current state of the solution and the solver to perform this evaluation. 

The function first updates the gradient and Hessian approximation. Then, it computes the predicted decrease in the objective function value (`Delta_mk`) and the actual decrease (`Delta_fk`). Based on these values, it calculates the quality ratio `ρk`.

# Errors
- Throws an error if `Delta_mk` is negative, which suggests an issue with the model's prediction.
- Throws an error if `Delta_mk` is zero, as it indicates no predicted improvement from the step, which is problematic for calculating the quality index value.
- A general error labeled `@error "floc"` is included as a catch-all for any unexpected scenarios.

"""
function checkQualityOfCandidateStep(
    SR1params, 
    pr,
    solState::SolStateType,
    pj::Vector{Float64},
    solverState::SolverStateType)

    @unpack xk, fk = solState
    @unpack gkm1, Bkm1 = SR1params
    gk, Bk = gkm1, Bkm1
    mk(pk) = fk + gk'*pk + 1/2*pk'*Bk*pk

    @unpack fevals = solverState

    f = pr.objective
    p = pr.p
    Delta_mk = fk - mk(pj)
    fj = f(xk+pj, p, getGradientToo=false)
    Delta_fk = fk - fj; fevals += 1

    @pack! solverState = fevals

    if Delta_mk < 0
        @error "problematic prediction which shouldn't even have been possible!"
    elseif Delta_mk == 0
        @error "zero Delta_mk? I'll never get my quality index value."
    else
        ρk = Delta_fk/Delta_mk
        return (ρk=ρk, fj=fj, solverState=solverState)
    end
    
    @error "floc"
end

