mutable struct SolStateType
    k::Int
    fk::Float64
    gk::Vector{Float64}
    gkmag::Float64
    pk::Vector{Float64}
    Hk::Matrix{Float64}
    alphak::Float64
    SolStateType() = new(0, 0.0, Float64[], 0.0, Float64[], Matrix{Float64}(undef, 0, 0), 0.0)

end

mutable struct SolverStateType
    k::Int
    fevals::Int
    gevals::Int
    Hevals::Int
    alpha_evals::Int # only current evals
    success_ls::Bool
    SolverStateType() = new(0, 0, 0, 0, 0, false)
end

# using Parameters
# solState = SolStateType()
# @unpack fk, gk, Hk = solState
# println("$fk, $gk, $Hk")
# fk = 3.2; 
# gk = [-1.0, 1.0]; 
# Hk = Float64.([5 -5;-5 5]);
# @pack! solState = fk, gk, Hk
# @unpack fk, gk, Hk = solState
# println("$fk, $gk, $Hk")

# solverState = SolverStateType()