include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("objfuns/functionQP.jl")
include("optimizeECQP.jl")
include("solveLP.jl")

function getECQPStep(prASQP, solStateASQP;
    verbose::Bool=false,
    verbose_ls::Bool=false,
)

    @unpack etol, itol = pr.alg
    @unpack objective, p, objectiveString = prASQP
    objectiveStringECQP = objectiveString*"_ECQPsubroutine"
    problemTypeECQP = "ECQP"
    methodECQP = "ProjectedGradientCG"

    pDictASQP = p
    @unpack G, c, mE, Ae, be, A, b = pDictASQP
    @unpack xk, Wk = solStateASQP

    myprintln(verbose, "We start for a good feasible descent direction from current point xk = $(xk)")
    WIk = Wk[mE+1:end]
    Ae = vcat(Ae, A[WIk .- mE, :])
    be = vcat(be, b[WIk .- mE])

    myprintln(verbose, "Let's look at the Ae and be being inserted into ECQP solver in order to represent our current Wk:")
    @show Ae, be

    myprintln(verbose, "Quick safety check to make sure we're inserting a feasible xk into our ECQP: What's Ae*xk - be?")
    rk = Ae*xk - be
    if norm(rk) < etol
        myprintln(verbose, "It's a zero vector as it should be.")
    else
        error("The residual is non-zero!. Then ECQP will give infeasible results as well.")
    end

    subroutineCall = true
    # @show subroutineCall
    pDictECQP = Dict(:G=>G, :c=>c, :Ae=>Ae, :be=>be, :subroutineCall=>subroutineCall)
    # pDictECQP = @packDict "{G, c, Ae, be}"
    # pDictECQP[:subroutineCall] = subroutineCall
    @show pDictECQP
    # @show pDictECQP
    prECQP = generate_pr(objective, xk, problemType=problemTypeECQP, method=methodECQP, params=pDictECQP, objectiveString=objectiveStringECQP, verbose=false)

    @show prECQP.p
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

function computeLagrangianMultipliersQP(xk, G, c, Awk)
    rhs = G*xk + c
    lambdask = transpose(Awk) \ rhs
    return lambdask
end
# c = [1, 3, 5, 2]
# n = length(c)
# x0 = rand(n)
# Ae = [1 1 9 5; 3 5 0 8; 2 0 6 13]
# be = [7, 3, 5]
# mE = length(be)
# A = [0 0 -1 -1]
# b = [-1]
# mI = length(b)
# pLP = Dict(:c => c, :Ae => Ae, :be => be, :mE => mE, :A => A, :b => b, :mI => mI)
# w0 = x0

# objective = LPObjectiveFunction
# # objectiveOriginal = LPObjectiveFunction
# # objectiveString = string(objectiveOriginal)
# objectiveString = "lpTestFunction1"
# params = pLP

# pr = generate_pr(objective, w0, params=params, problemType="LP"; objectiveString=objectiveString)

# xf = computeFeasiblePointForLinearConstraints(pr.p)

# xopt, fopt = solveLP(pr)
