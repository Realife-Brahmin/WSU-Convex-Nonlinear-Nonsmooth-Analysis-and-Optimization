include("objective.jl")
include("../solveLP.jl")

using HiGHS
using JuMP
using Parameters

function LPObjectiveFunction(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack c = pDict
    f = transpose(x)*c

    if !getGradientToo
        return f
    elseif getGradientToo
        g = c
        return f, g
    else
        @error "floc"
    end
    
end

# c = [1, 3, 5, 2]
# n = length(c)
# x0 = rand(n)
# Ae = [1 1 9 5; 3 5 0 8; 2 0 6 13]
# be = [7, 3, 5]
# A = [0 0 -1 -1]
# b = [-1]
# pLP = Dict(:c => c, :Ae => Ae, :be => be, :A => A, :b => b)
# w0 = x0
# vector_model = Model(HiGHS.Optimizer)
# @variable(vector_model, x[1:n] >= 0)
# @constraint(vector_model, Ae * x .== be)
# @constraint(vector_model, A * x .>= b)
# @objective(vector_model, Min, c' * x)
# optimize!(vector_model)
# @assert is_solved_and_feasible(vector_model)
# xopt = value.(x)
# fopt = objective_value(vector_model)
# objective = LPObjectiveFunction
# # objectiveOriginal = LPObjectiveFunction
# # objectiveString = string(objectiveOriginal)
# objectiveString = "lpTestFunction1"
# params = pLP

# pr = generate_pr(objective, w0, params=params, problemType="LP"; objectiveString=objectiveString)

# xopt, fopt = solveLP(pr)
# f0 = equalityConstrainedQP(w0, pECQP)

# solState = SolStatePGCGType(w0, G, c, A, fk=f0)
