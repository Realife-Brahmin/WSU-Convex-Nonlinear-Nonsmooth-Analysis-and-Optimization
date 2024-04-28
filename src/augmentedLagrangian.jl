include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")
include("optimize2.jl") # to call QN-BFGS

using LinearAlgebra
using Parameters

function solveAugmentedLagrangianFunction(prALP, solStateALP;
    verbose::Bool=false,
    verbose_ls::Bool=false,
)

    @unpack etol, itol = pr.alg
    @unpack objective, p, objectiveString = prALP
    objectiveStringECQP = objectiveString * "_ECQPsubroutine"
    problemTypeECQP = "ECQP"
    methodECQP = "ProjectedGradientCG"

    pDictALP = p
    @unpack G, c, mE, Ae, be, A, b = pDictALP
    @unpack xk, Wk = solStateALP

    myprintln(verbose, "We start for a good feasible descent direction from current point xk = $(xk)")
    WIk = Wk[mE+1:end]
    Ae = vcat(Ae, A[WIk.-mE, :])
    be = vcat(be, b[WIk.-mE])

    # myprintln(verbose, "Let's look at the Ae and be being inserted into ECQP solver in order to represent our current Wk:")
    # @show Ae, be

    subroutineCall = true
    # @show subroutineCall
    pDictECQP = Dict(:G => G, :c => c, :Ae => Ae, :be => be, :subroutineCall => subroutineCall)
    # pDictECQP = @packDict "{G, c, Ae, be}"
    # pDictECQP[:subroutineCall] = subroutineCall
    # @show pDictECQP
    # @show pDictECQP
    prECQP = generate_pr(objective, xk, problemType=problemTypeECQP, method=methodECQP, params=pDictECQP, objectiveString=objectiveStringECQP, verbose=false)

    # @show prECQP.p
    res = optimizeECQP(prECQP, verbose=verbose, verbose_ls=verbose_ls, log=false)

    xvals = res[:xvals]
    itr = size(xvals, 2)
    if itr == 0 # x0 is xkp1
        xkp1 = xk
    else
        xkp1 = xvals[:, itr]
    end

    pk = xkp1 - xk
    return pk

end

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

