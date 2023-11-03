using Parameters

include("helperFunctions.jl")
include("types.jl")

function StrongWolfe(pr::NamedTuple, 
    solState::SolStateType, 
    solverState::SolverStateType;

    alphaLo::Float64=0.0, alphaHi::Float64=100.0, 
    alphatol::Float64=1e-10,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"

    function myprintln1(v, message)
        myprintln(v, message, log=log, log_path=log_txt)
    end

    c1 = pr.alg.c1
    c2 = pr.alg.c2
    obj = pr.objective
    p = pr.p

    @unpack xkm1, xk, fkm1 ,fk, gkm1, gk, pkm1, pk = solState
    
    fevals = 0; gevals = 0
    interpolParams = InterpolParams(j=1, alphaj=alphaHi, alphaHi=alphaHi, alphaLo=alphaLo, alphatol=alphatol)

    keepSearching = true
    while keepSearching

        @unpack alphaj = interpolParams

        xj = xk + alphaj*pk
        fj = obj(xj, p, getGradientToo=false)
        fevals += 1

        if StrongWolfe1(fk, fj, gk, pk, alphaj)
            myprintln1(true, )
        end

    end

    return (solState=solState, solverState=solverState)
end

function StrongWolfe1(fk, fj, gk, pk, alphaj; c1=1e-4)
    if fk - fj ≥ c1*alphaj*gk'*pk
        return true
    else
        return false
    end

    @error "floc"
end 

function StrongWolfe2(gk, gj, pk, alphaj; c2=0.9)
    if abs(gj'*pk) ≤ c2*abs(gk'*pk)
        return true
    else
        return false
    end

    @error "floc"
end