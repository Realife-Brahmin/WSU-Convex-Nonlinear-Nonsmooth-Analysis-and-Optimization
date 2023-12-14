using LinearAlgebra
using Parameters
using Test

include("helperFunctions.jl")
include("types.jl")

"""
    checkQualityOfCandidateStep(SR1params, pr, xk, fk, pj, solverState)

Evaluate the quality of a candidate step in a numerical optimization algorithm.

# Arguments
- `SR1params`: Parameters for SR1 (Symmetric Rank 1) update, often used in quasi-Newton methods.
    It should contain at least `gk` (gradient at `xk`) and `Bk` (approximation to the Hessian matrix).
- `pr`: Problem-related structure, specifics of which depend on the problem being solved.
- `xk::Vector{Float64}`: Current point in the optimization space.
- `fk::Float64`: Function value at `xk`.
- `pj::Vector{Float64}`: The step to be evaluated, proposed from `xk`.
- `solverState::SolverStateType`: Contains state information of the solver, like function evaluation count.

# Returns
- A tuple with the ratio `ρk` (indicating the quality of the step) and the updated `solverState`.

# Description
This function assesses the quality of a candidate step (`pj`) in an optimization process. 
It uses the SR1 update parameters and the current state of the solver to evaluate the step.
If the step improves the solution, it calculates the ratio `ρk` of actual decrease in the 
objective function to the predicted decrease. If `ρk` is significantly positive, it implies 
a good quality step. The function also updates the solver state, particularly the count of 
function evaluations.

# Errors
Throws an error if the predicted decrease (`Delta_mk`) is negative or zero, indicating 
anomalies in the optimization process.
"""
function checkQualityOfCandidateStep(SR1params, pr,
    xk::Vector{Float64},
    fk::Float64,
    pj::Vector{Float64},
    solverState::SolverStateType)

    @unpack gkm1, Bkm1 = SR1params
    gk, Bk = gkm1, Bk
    mk(pk) = fk + gk'*pk + 1/2*pk'*Bk*pk

    @unpack fevals = solverState

    f = pr.obj
    p = pr.p
    Delta_mk = fk - mk(pj)
    Delta_fk = fk - f(xk+pj, p, getGradientToo=false); fevals += 1

    @pack! fevals = solverState

    if Delta_mk < 0
        @error "problematic prediction which shouldn't even have been possible!"
    elseif Delta_mk == 0
        @error "zero Delta_mk? I'll never get my quality index value."
    else
        ρk = Delta_fk/Delta_mk
        return (ρk=ρk, solverState=solverState)
    end
    
    @error "floc"
end

