using HiGHS
using JuMP
using Parameters

function solveLP(pr;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    @unpack c, Ae, be, A, b = pr.p[:params]
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