using Parameters

include("helperFunctions.jl")
include("interpolation.jl")
include("types.jl")

function StrongWolfe(pr::NamedTuple, 
    solState::SolStateType, 
    solverState::SolverStateType;

    alphaLo::Float64=0.0, 
    alphaHi::Float64=1.0,
    alpha0=nothing,
    alphatol::Float64=1e-13,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg[:method]*"_"*pr.alg[:linesearch]*"_"*string(pr.alg[:maxiter])*".txt"

    function myprintln1(v, message)
        myprintln(v, message, log=log, log_path=log_txt)
    end

    c1 = pr.alg[:c1]
    c2 = pr.alg[:c2]
    obj = pr.objective
    p = pr.p

    if isnothing(alpha0)
        alphaj = min((alphaLo + alphaHi) / 2, 1.0)
    else
        alphaj = alpha0
    end

    @unpack xkm1, xk, fkm1 ,fk, gkm1, gk, pkm1, pk = solState
    
    @unpack fevals, gevals = solverState
    interpolParams = InterpolParams(j=1, alphaj=alphaj, alphaHi=alphaHi, alphaLo=alphaLo, alphatol=alphatol)

    keepSearching = true # for the while loop
    success_ls = false # a field of solverState useful outside
    SW1_satisfied_once = false # to cover a special case where the initial point is already optimal
    

    j = 1
    xj = xk; fj = fk; gj = gk
    # @show xj, fj, fk, gj
    while keepSearching

        @unpack alphaj = interpolParams
        xj = xk + alphaj*pk
        fj = obj(xj, p, getGradientToo=false)
        fevals += 1

        # verbose = true

        if StrongWolfe1(fk, fj, gk, pk, alphaj)
            myprintln1(verbose, "SW1 satisfied for αⱼ = $(alphaj)")
            SW1_satisfied_once = true
            @show fj, gj = obj(xj, p, getGradientToo=true)
            fevals += 1; gevals += 1
            if StrongWolfe2(gk, gj, pk)
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
            @warn "LineSearch failed. Returning best obtained xⱼ obtained at αⱼ = $(alphaj)"
            keepSearching = false
        end

        j += 1
    end

    @unpack alphaj = interpolParams
    xkm1 = xk
    fkm1 = fk
    gkm1 = gk
    pkm1 = pk
    gmagkm1 = sum(abs.(gk))
    @pack! solState = xkm1, fkm1, gkm1, pkm1, gmagkm1
    if SW1_satisfied_once
        alphak = alphaj
        xk = xj
        fk = fj
        gk = gj
        gmagk = sum(abs.(gj))
        @pack! solState = xk, fk, gk, gmagk, alphak
    else
        @warn "Function didn't decrease for any step length. Is it already optimal?"
    end

    alpha_evals = j
    @pack! solverState = fevals, gevals, alpha_evals, success_ls

    return (solState=solState, solverState=solverState)
end

function StrongWolfe1(fk, fj, gk, pk, alphaj; c1=1e-4)
    if fj - fk < c1*alphaj*gk'*pk
        return true
    else
        return false
    end

    @error "floc"
end 

function StrongWolfe2(gk, gj, pk, c2=0.9)
    if gj'*pk ≥ c2*gk'*pk
        return true
    else
        return false
    end

    @error "floc"
end