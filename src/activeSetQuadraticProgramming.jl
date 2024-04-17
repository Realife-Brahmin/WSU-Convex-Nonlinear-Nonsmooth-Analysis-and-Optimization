include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("solveLP.jl")

using HiGHS
using JuMP
using Parameters

function computeFeasiblePointForLinearConstraints(pDict;
    )

    @unpack mE, mI, Ae, be, A, b = pDict[:params]

    n = size(Ae, 2)

    vector_model = Model(HiGHS.Optimizer)
    @variable(vector_model, x[1:n])
    @constraint(vector_model, Ae * x .== be)
    @constraint(vector_model, A * x .>= b)
    optimize!(vector_model)
    @assert is_solved_and_feasible(vector_model)
    xopt = value.(x)

    return xopt
    
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
