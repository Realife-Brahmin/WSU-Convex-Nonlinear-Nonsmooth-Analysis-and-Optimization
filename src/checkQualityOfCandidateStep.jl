using LinearAlgebra
using Parameters
using Test

include("helperFunctions.jl")
include("types.jl")

function checkQualityOfCandidateStep(mk::Function,
    f::Function,
    xk::Vector{Float64},
    fk::Float64,
    pj::Vector{Float64},
    solverState::SolverStateType)

    @unpack fevals = solverState

    Delta_mk = fk - mk(pj)
    Delta_fk = fk - f(xk+pj); fevals += 1

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

