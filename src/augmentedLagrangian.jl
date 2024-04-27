include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")
include("optimize2.jl") # to call QN-BFGS

function ALOBJ(w, p;
    getGradientToo::Bool=false)

    n = length(w) - p.mI
    x = w[1:n]
    y = w[n+1:end]

    f, g = p.obj(x, p)
    cE, hE = p.econ(x, p)
    cI, hI = p.icon(x, p)

    Y = diag(y)
    c = vcat(cE, cI-Y*Y)
    F = f - transpose(lambdas)*c + mu/2 * transpose(c) * c

    if !getGradientToo
        return F
    elseif getGradientToo
        G = vcat(g, zeros(p.mI)) - vcat( hcat(cE, cI), zeros(p.mI, n) - 2*Y) * (lambdas - mu*c)
        return F, G
    else
        @error("floc")
    end

    @error("floc")
end