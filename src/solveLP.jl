using HiGHS
using JuMP
using Parameters

include("helperFunctions.jl")


function computeFeasiblePointForLinearConstraints(pDict;
    )

    @unpack lb, ub, mE, mI, Ae, be, A, b = pDict

    n = size(Ae, 2)
    isempty(lb) ? lb = zeros(n) : lb = lb
    isempty(ub) ? ub = myfill(lb, Inf) : ub = ub

    model = Model(HiGHS.Optimizer)
    @variable(model, lb[i] .<= x[i=1:n] .<= ub[i])
    @constraint(model, Ae * x .== be)
    @constraint(model, A * x .>= b)
    optimize!(model)
    @assert is_solved_and_feasible(model)
    xfeas = value.(x)

    return xfeas
    
end

"""
    solveLP(pDict; verbose=false, verbose_ls=false, log=true, log_path="./logging/")

Solve a linear programming (LP) problem using the HiGHS optimizer.

# Arguments
- `pDict::Dict`: A dictionary containing the parameters of the LP problem. Expected keys are:
    - `c`: Coefficient vector for the objective function.
    - `lb`: Lower bounds for decision variables. If not provided, default to zero.
    - `ub`: Upper bounds for decision variables. If not provided, default to `Inf`.
    - `Ae`: Coefficient matrix for equality constraints.
    - `be`: Right-hand side vector for equality constraints.
    - `A`: Coefficient matrix for inequality constraints.
    - `b`: Right-hand side vector for inequality constraints.

# Keyword Arguments
- `verbose::Bool=false`: Flag to enable verbose output.
- `verbose_ls::Bool=false`: Flag to enable verbose output for line search (not used in current implementation).
- `log::Bool=true`: Flag to enable logging.
- `log_path::String="./logging/"`: Path to the directory where logs should be saved.

# Returns
- `xopt`: The optimal solution vector.
- `fopt`: The optimal value of the objective function.

# Notes
- `is_solved_and_feasible` is an assert statement that checks if the model has been solved and is feasible.
- The function uses the HiGHS solver, which must be installed and included in the project environment.
- Box constraints are applied to the decision variables based on provided `lb` and `ub` vectors.
- The model allows decision variables to take on any real value within the specified bounds.

# Example
```julia
# Define problem parameters
c = [1, 3, 5, 2]
lb = [0, 0, 0, 0]  # Lower bounds for variables
ub = [Inf, Inf, Inf, Inf]  # Upper bounds for variables
Ae = [1 1 9 5; 3 5 0 8; 2 0 6 13]
be = [7, 3, 5]
A = [0 0 -1 -1]
b = [-1]

# Pack LP problem parameters into a dictionary
params = @packDict "{c, lb, ub, Ae, be, A, b}"
pLP = (params=params, data=[])

# Solve the LP problem
xopt, fopt = solveLP(pLP)
println("Optimal solution: ", xopt)
println("Optimal objective value: ", fopt)
```
"""
function solveLP(pDict;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    @show pDict[:params]
    @unpack c, lb, ub, Ae, be, A, b = pDict[:params]
    n = length(c)

    isempty(lb) ? lb = zeros(n) : lb = lb
    isempty(ub) ? ub = myfill(lb, Inf) : ub = ub

    vector_model = Model(HiGHS.Optimizer)
    @variable(vector_model, lb[i] .<= x[i=1:n] .<= ub[i])  # Define the variables with box constraints
    @constraint(vector_model, Ae * x .== be)
    @constraint(vector_model, A * x .>= b)
    @objective(vector_model, Min, c' * x)
    optimize!(vector_model)
    @assert is_solved_and_feasible(vector_model)
    xopt = value.(x)
    fopt = objective_value(vector_model)

    return xopt, fopt

end

# c = [1, 3, 5, 2]
# n = length(c)
# x0 = rand(n)
# lb = myzeros(x0)
# ub = myfill(x0, Inf)
# Ae = [1 1 9 5; 3 5 0 8; 2 0 6 13]
# be = [7, 3, 5]
# mE = length(be)
# A = [0 0 -1 -1]
# b = [-1]
# mI = length(b)

# params = @packDict "{c, lb, ub, mE, Ae, be, mI, A, b}"
# pLP = (params=params, data=[])
# xopt, fopt = solveLP(pLP)
