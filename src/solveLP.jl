using HiGHS
using JuMP
using Parameters

include("helperFunctions.jl")

function solveLP(pDict;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    @show pDict[:params]
    @unpack c, Ae, be, A, b = pDict[:params]
    n = length(c)

    vector_model = Model(HiGHS.Optimizer)
    @variable(vector_model, x[1:n] >= 0)
    @constraint(vector_model, Ae * x .== be)
    @constraint(vector_model, A * x .>= b)
    @objective(vector_model, Min, c' * x)
    optimize!(vector_model)
    @assert is_solved_and_feasible(vector_model)
    xopt = value.(x)
    fopt = objective_value(vector_model)

    return xopt, fopt

end

c = [1, 3, 5, 2]
n = length(c)
x0 = rand(n)
Ae = [1 1 9 5; 3 5 0 8; 2 0 6 13]
be = [7, 3, 5]
mE = length(be)
A = [0 0 -1 -1]
b = [-1]
mI = length(b)

params = @packDict "{c, mE, Ae, be, mI, A, b}"
pLP = (params=params, data=[])
xopt, fopt = solveLP(pLP)
