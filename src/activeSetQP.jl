include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("objfuns/functionQP.jl")
include("projectedGradientConjugateGradient.jl")
include("solveLP.jl")

function getECQPStep(prASQP, pDict;
    verbose::Bool=false,
)
    error("Okay I still need to write a subroutine format of optimizeECQP")
    solverState = SolverStateECQPType()

    # @unpack genetic parameters = pr.alg

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    etol = pr.alg[:etol]

    return pk
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
