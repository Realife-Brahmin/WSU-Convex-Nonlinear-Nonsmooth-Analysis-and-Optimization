using Parameters

include("helperFunctions.jl")
include("interpolation.jl")
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
    
    @unpack fevals, gevals = solverState
    interpolParams = InterpolParams(j=1, alphaj=alphaHi, alphaHi=alphaHi, alphaLo=alphaLo, alphatol=alphatol)

    keepSearching = true # for the while loop
    success_ls = false # a field of solverState useful outside

    j = 1
    xj = xk; fj = fk; gj = gk
    
    while keepSearching

        @unpack alphaj = interpolParams
        @show xk, alphaj, pk
        xj = xk + alphaj*pk
        fj = obj(xj, p, getGradientToo=false)
        fevals += 1

        if StrongWolfe1(fk, fj, gk, pk, alphaj)
            myprintln1(verbose, "SW1 satisfied for αⱼ = $(alphaj)")
            fj, gj = obj(xj, p)
            fevals += 1; gevals += 1
            if StrongWolfe2(gk, gj, pk, alphaj)
                myprintln1(verbose, "SW2 satisfied for αⱼ = $(alphaj)")
                success_ls = true
                keepSearching = false
            else
                myprintln1(verbose, "SW2 NOT satisfied. αⱼ ↑")
                interpolParams.change = "increase"
                interpolParams = bisection(interpolParams)
            end
        else
            myprintln(verbose, "SW1 NOT satisfied. αⱼ ↓")
            interpolParams.change = "decrease"
            interpolParams = bisection(interpolParams)
        end

        @unpack alphatolBreached = interpolParams
        if alphatolBreached
            @warn "LineSearch failed. Returning best obtained xⱼ."
            keepSearching = false
        end

        j += 1
    end

    @unpack alphaj = interpolParams
    xkm1 = xk
    fkm1 = fk
    gkm1 = gk
    pkm1 = pk
    @pack! solState = xkm1, fkm1, gkm1, pkm1
    alphak = alphaj
    xk = xj
    fk = fj
    gk = gj
    gmagk = sum(abs.(gj))
    @pack! solState = xk, fk, gk, gmagk, alphak
    alpha_evals = j
    @pack! solverState = fevals, gevals, alpha_evals, success_ls

    println(solState)
    println(solverState)
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