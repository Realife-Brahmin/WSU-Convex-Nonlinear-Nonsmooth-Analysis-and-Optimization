include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")
include("optimize2.jl") # to call QN-BFGS

using LinearAlgebra
using Parameters

function ALOBJ(w, psubDict;
    getGradientToo::Bool=false)

    @unpack lambda, mu, mE, mI = psubDict
    n = length(w) - mI # length of x
    x = w[1:n]
    y = w[n+1:end]

    f, g = psubDict.objective(x, psubDict, getGradientToo=true)
    cE = psubDict.econ(x, psubDict, getGradientToo=false)
    cI = psubDict.icon(x, psubDict, getGradientToo=false)

    # Y = diagm(y)
    y2 = transpose(y)*y
    c = vcat(cE, cI - y2)
    F = f - transpose(lambda)*c + mu/2 * transpose(c) * c

    if !getGradientToo
        return F
    elseif getGradientToo
        G = vcat(g, zeros(mI)) - vcat( hcat(cE, cI), zeros(mI, n) - 2*y2) * (lambda - mu*c) # where are lambda and mu coming from
        return F, G
    else
        @error("floc")
    end

    @error("floc")
end