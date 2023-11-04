using Parameters

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
    

    function SolStateType(;
        k=1, 
        xkm1=Float64[], xk=Float64[],
        fkm1=0.0, fk=0.0,
        gkm1=Float64[], gk=Float64[],
        gmagkm1=0.0, gmagk=0.0, 
        pkm1=Float64[], pk=Float64[], 
        Hkm1=Matrix{Float64}(undef, 0, 0), Hk=Matrix{Float64}(undef, 0, 0), 
        alphak=0.0)

        new(k, xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk, pkm1, pk, Hkm1, Hk, alphak)
    end

end

mutable struct SolverStateType
    k::Int
    fevals::Int
    gevals::Int
    Hevals::Int
    alpha_evals::Int # only current evals
    success_ls::Bool

    function SolverStateType(;
        k=1, 
        fevals=0, 
        gevals=0, 
        Hevals=0, 
        success_ls=false)
        new(k, fevals, gevals, Hevals, success_ls)
    end
end

mutable struct InterpolParams
    j::Int
    alphaj::Float64
    alphaLo::Float64
    alphaHi::Float64
    alphatol::Float64
    alphatolBreached::Bool
    change::String

    function InterpolParams(;j=1, alphaj=100.0, alphaLo=0.0, alphaHi=100.0, alphatol=1e-10, alphatolBreached=false, change="noChange")
        new(j, alphaj, alphaLo, alphaHi, alphatol, alphatolBreached, change)
    end

end

# solState = SolStateType()
# solverState = SolverStateType()
# interpolParams = InterpolParams(alphatol=33)