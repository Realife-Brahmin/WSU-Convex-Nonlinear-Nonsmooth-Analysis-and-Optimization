include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("objfuns/functionQP.jl")
include("optimizeECQP.jl")
include("solveLP.jl")

function getECQPStep(prASQP, solStateASQP;
    verbose::Bool=false,
    verbose_ls::Bool=false,
)

    @unpack objective, p, objectiveString = prASQP
    objectiveStringECQP = objectiveString*"_ECQPsubroutine"
    problemTypeECQP = "ECQP"
    methodECQP = "ProjectedGradientCG"

    pASQP = p
    @unpack G, c, mE, Ae, be, A, b = pASQP
    @unpack xk, Wk = solStateASQP

    WIk = Wk[mE+1:end]
    Ae = vcat(Ae, A[WIk .- mE, :])
    be = vcat(be, b[WIk .- mE])

    subroutineCall = true
    # @show subroutineCall
    pECQP = @packDict "{G, c, Ae, be}"
    pECQP[:subroutineCall]=subroutineCall

    # @show pECQP
    prECQP = generate_pr(objective, xk, problemType=problemTypeECQP, method=methodECQP, params=pECQP, objectiveString=objectiveStringECQP, verbose=false)
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
