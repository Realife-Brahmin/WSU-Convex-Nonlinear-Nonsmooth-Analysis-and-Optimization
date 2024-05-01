include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")
include("optimize2.jl") # to call QN-BFGS

using LinearAlgebra
using Parameters

function combinedConstraintsALP(wk, pDictUnc)
    @unpack n, mE, mI, econ, icon = pDictUnc
    xk, yk = wk[1:n], wk[n+1:end]
    cE = econ(xk, pDictUnc, getGradientToo=false)
    cI = icon(xk, pDictUnc, getGradientToo=false)
    c = vcat(cE, cI - transpose(yk)*yk)
    return c
end

function solveAugmentedLagrangianFunction(prALP, solStateALP;
    verbose::Bool=false,
    verbose_ls::Bool=false,
)

    @unpack etol, mutol, dxtol, gtol = pr.alg
    @unpack objective, p, objectiveString = prALP
    objectiveALP = objective
    objectiveUnc = ALOBJ
    objectiveStringUnc = objectiveString * "_ALPsubroutine"
    problemTypeUnc = "Unconstrained"
    methodUnc = "QuasiNewton"

    pDictALP = p # p is indeed pDict from the original problem
    @unpack mE, econ, mI, icon  = pDictALP
    @unpack wk, lambdak, muk, tk = solStateALP # again, wk is actually wk

    myprintln(verbose, "From the current point wk = $(wk)")

    subroutineCall = true

    addendum = Dict(:subroutineCall => subroutineCall, :lambda => lambdak, :mu => muk, :objective => objectiveUnc, :objectiveALP => objectiveALP, :tau => tk)

    pDictUnc = merge(deepcopy(pDictALP), addendum)

    # @show pDictUnc

    prUnc = generate_pr(objectiveUnc, wk, problemType=problemTypeUnc, method=methodUnc, params=pDictUnc, objectiveString=objectiveStringUnc, verbose=false)

    # @show prUnc

    # error("Stop here.")
    resUnc = optimize2(prUnc, verbose=verbose, verbose_ls=verbose_ls, log=false)

    @unpack xopt, fopt, iter = resUnc
    coptUnc = combinedConstraintsALP(xopt, pDictUnc) #do # constraint violations
    return xopt, fopt, coptUnc, iter

end

function ALOBJ(w, psubDict;
    getGradientToo::Bool=false)

    @unpack n, m, lambda, mu, mE, mI, econ, icon, objectiveALP = psubDict
    # @show n, w 
    # n = length(w) - mI # length of x
    x = w[1:n]
    y = w[n+1:end]
    # error("Okay we're done.")

    @show x, psubDict
    
    f, g = objectiveALP(x, psubDict, getGradientToo=true)
    cE = econ(x, psubDict, getGradientToo=false)
    cI = icon(x, psubDict, getGradientToo=false)

    c = vcat(cE, cI - y.^2) # mE+mI length c vector
    F = f - transpose(lambda)*c + 1//2 * mu * transpose(c) * c

    if !getGradientToo
        return F
    elseif getGradientToo
        _ , hE = econ(x, psubDict, getGradientToo=true)
        _ , hI = icon(x, psubDict, getGradientToo=true)
        h = vcat(hE, hI) # mxn matrix
        Y = diagm(y) # mI x mI diagonal matrix to be used in gradient of inequality constraints
        hy = vcat( zeros(mE, mI), -2 * Y)
        dL1dw = vcat( g, zeros(mI) ) # (Nx1) vector # derivative of first lagrangian term (f(x) wrt w)
        dL23dw = vcat(transpose(h), transpose(hy)) * (-lambda + mu * c) # (Nx1) vector # concatenated derivatives of second and third lagrangian terms (-λ*c(w) + 1//2*μ*transpose(c(w))*c(w)) wrt w)
        G = dL1dw + dL23dw
        return F, G
    else
        @error("floc")
    end

    @error("floc")
end

