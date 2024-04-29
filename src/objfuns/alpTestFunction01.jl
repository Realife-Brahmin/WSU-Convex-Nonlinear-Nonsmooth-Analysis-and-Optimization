include("objective.jl")
include("../AugmentedLagrangian.jl")

using Parameters

function cE01(x, p; getGradientToo::Bool=false)
    mE = 1
    cE = zeros(mE)
    cE[1] = x[1]^2 + x[2]^2 - 1

    if !getGradientToo
        return cE
    elseif getGradientToo
        hE = [2*x[1], 2*x[2]]
        return cE, hE
    else
        @error("floc")
    end

    @error("floc")

end

mE = 1

function cI01(x, p; getGradientToo::Bool=false)
    mI = 2
    cI = zeros(mI)
    cI[1] = x[1]^3 + x[2]
    cI[2] = x[1]^2 + 2*x[2]^2 - 3
    if !getGradientToo
        return cI
    elseif getGradientToo
        n = length(x)
        hI = zeros(mI, n)
        hI[1, 1] = 3*x[1]^2
        hI[1, 2] = 1
        hI[2, 1] = 2*x[1]
        hI[2, 2] = 4*x[2]
        return cI, hI
    else
        @error("floc")
    end

    @error("floc")

end

mI = 2

function alpTestFunction01(x, p, getGradientToo::Bool=true)
    f = (x[1]-1)^2 + 100*(x[2]-x[1]^2)^2
    if !getGradientToo
        return f
    elseif getGradientToo
        g = [2*(x[1]-1) - 400*x[1]*(x[2]-x[1]^2), 200*(x[2]-x[1]^2)]
        return f, g
    else
        @error("floc")
    end

    @error("floc")
end

objective = alpTestFunction01
objectiveOriginal = alpTestFunction01
objectiveString = "alpTestFunction01"
problemType = "Constrained"
pALP = Dict(:mE=>mE, :econ=>cE01, :mI=>mI, :icon=>cI01)
# pr = generate_pr(objective, x0, params=pALP, problemType=problemType; objectiveString=objectiveString)

using JuMP, Ipopt

model = Model(Ipopt.Optimizer)
x0 = [-1.2, 1.0]
n = length(x0)
y0 = zeros(mI)
# @variable(model, w[i=1:4], start = [1.21, 1.79, 1.1, sqrt(0.79)][i])
@variable(model, x[i=1:n], start = x0[i])
@variable(model, y[i=1:mI], start = y0[i])

@objective(model, Min, (x[1] - 1)^2 + (x[2] - 2)^2)

# @objective(model, Min, (w[1] - 1)^2 + (w[2] - 2)^2)

# Constraint equations assuming w[3] and w[4] are slack variables for the inequalities
# w[1]^2 <= 1 => w[1]^2 - w[3]^2 = 0 (slack variable w[3] squared)
# w[2] - 1 <= 0 => w[2] - 1 - w[4]^2 = 0 (slack variable w[4] squared)
# We also have an equality constraint w[1] + w[2] = 5
# cI[1] = x[1]^3 + x[2]
# cI[2] = x[1]^2 + 2 * x[2]^2 - 3
# Rewrite the constraints to use the slack variables
# @constraint(model, w[1] - w[3]^2 == 0) # Equality constraint with slack variable w[3]
@constraint(model, x[1] + x[2] - 3 == 0)

# @constraint(model, x[1] - y[1]^2 == 0)
# @constraint(model, x[2] - 1 - y[2]^2 == 0)

# @constraint(model, x[1] + w[2] == 5) # Regular equality constraint
# @constraint(model, x[1]^2 + 2*x[2]^2 - y[2]^2 == 0) # Regular equality constraint

optimize!(model)

# w_optimal = [value(w[i]) for i in 1:4]
xopt = [value(x[i]) for i ∈ 1:n]
yopt = [value(y[i]) for i ∈ 1:mI]
f_optimal = objective_value(model)

# Call the function to solve the problem
# println("Optimal w: ", w_optimal)
println("Optimal Variables x: ", xopt)
println("Optimal Slack Variables y: ", yopt)
println("Optimal objective value: ", f_optimal)


# optimize!(nlp)
# println(solution_summary(nlp))
# wopt = [value(w[i]) for i ∈ eachindex(w)]
# fopt = objective_value(nlp)
# myprintln(true, "Optimal decision variables are w = $(wopt)")
# myprintln(true, "Optimal objective value is f = $(fopt)")

