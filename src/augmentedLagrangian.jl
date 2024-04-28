include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")
include("optimize2.jl") # to call QN-BFGS

using LinearAlgebra
using Parameters

function ALOBJ(w, psubDict;
    getGradientToo::Bool=false)

    @unpack lambda, mu, mE, mI, econ, icon, objectiveALP = psubDict
    n = length(w) - mI # length of x
    x = w[1:n]
    y = w[n+1:end]

    f, g = objectiveALP(x, psubDict, getGradientToo=true)
    cE = econ(x, psubDict, getGradientToo=false)
    cI = icon(x, psubDict, getGradientToo=false)

    Y = diagm(y)
    # y2 = transpose(y)*y
    c = vcat(cE, cI - transpose(y) * y) # mE+mI vector
    F = f - transpose(lambda)*c + 1//2 * mu * transpose(c) * c

    if !getGradientToo
        return F
    elseif getGradientToo
        G = vcat(g, zeros(mI)) - vcat( hcat(cE, cI), zeros(mI, n) - 2*Y) * (lambda - mu*c) # where are lambda and mu coming from
        return F, G
    else
        @error("floc")
    end

    @error("floc")
end