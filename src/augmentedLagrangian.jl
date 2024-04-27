include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")
include("optimize2.jl") # to call QN-BFGS

using LinearAlgebra
using Parameters

function ALOBJ(w, psubDict;
    getGradientToo::Bool=false)

    @unpack lambdas, mu = psubDict
    n = length(w) - psubDict.mI # length of x
    x = w[1:n]
    y = w[n+1:end]

    f, g = psubDict.objective(x, psubDict, getGradientToo=true)
    cE = psubDict.econ(x, psubDict, getGradientToo=false)
    cI = psubDict.icon(x, psubDict, getGradientToo=false)

    Y = diagm(y)
    c = vcat(cE, cI-Y*Y)
    F = f - transpose(lambdas)*c + mu/2 * transpose(c) * c

    if !getGradientToo
        return F
    elseif getGradientToo
        G = vcat(g, zeros(psubDict.mI)) - vcat( hcat(cE, cI), zeros(psubDict.mI, n) - 2*Y) * (lambdas - mu*c) # where are lambdas and mu coming from
        return F, G
    else
        @error("floc")
    end

    @error("floc")
end