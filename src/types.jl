mutable struct SolStateType
    k::Int
    xkm1::Vector{Float64}
    xk::Vector{Float64}
    fkm1::Float64
    fk::Float64
    gkm1::Vector{Float64}
    gk::Vector{Float64}
    gmagkm1::Float64
    gmagk::Float64
    pkm1::Vector{Float64}
    pk::Vector{Float64}
    Hkm1::Matrix{Float64}
    Hk::Matrix{Float64}
    alphak::Float64
    
    SolStateType() = new(
        1, # k 
        Float64[], Float64[], # xkm1, xk
        0.0, 0.0,  # fkm1, fk
        Float64[], Float64[], # gkm1, gk 
        0.0, 0.0, # gmagkm1, gmagk
        Float64[], Float64[], # pkm1, pk  
        Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0), # Hkm1, Hk 
        0.0 # alphak
    )

end

mutable struct SolverStateType
    k::Int
    fevals::Int
    gevals::Int
    Hevals::Int
    alpha_evals::Int # only current evals
    success_ls::Bool

    SolverStateType() = new(1, 0, 0, 0, 0, false)
end

mutable struct InterpolParams
    j::Int
    alphaj::Float64
    alphaLo::Float64
    alphaHi::Float64
    alphatol::Float64
    alphatolBreached::Bool
    dir::String

    InterpolParams() = new(
        1, # j
        100, 0, 100, # alphaj, alphaLo, alphaHi
        1e-10, # alphatol
        false, # alphatolBreached
        "noChange" # dir
    )

end

# solState = SolStateType()
# solverState = SolverStateType()
interpolParams = InterpolParams()