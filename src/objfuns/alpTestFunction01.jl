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
# psubALP = deepcopy(pALP)
# lambda = zeros(mE+mI)
# mu = 1e1*ones(mE+mI)

# x0 = [-1.2, 1.0]
x0 = [3, 0]
n = length(x0)
# y0 = sqrt.(max.(0, cI01(x0, pALP, getGradientToo=false)))
y0 = zeros(2)
w0 = vcat(x0, y0)
using JuMP, NLPModelsJuMP, Percival
nlp = Model(NLPModelsJuMP.Optimizer)
set_attribute(nlp, "solver", Percival.PercivalSolver)
@variable(nlp, w[i=1:4], start = w0[i])
# @variable(nlp, x[i=1:n], start = x0[i])
# @variable(nlp, y[i=1:mI], start = y0[i])
@objective(nlp, Min, (w[1]-1)^2 + (w[2]-2)^2)
# @objective(nlp, Min, (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2)
# @objective(nlp, Min, (w[1] - 1)^2 + 100 * (w[2] - w[1]^2)^2)
# @constraint(nlp, w[1]^2 + w[2]^2 - 1 == 0)
@constraint(nlp, w[1] + w[2] - 3 == 0)
# @constraint(nlp, w[1]^3 + w[2] - w[3]^2 == 0 )
@constraint(nlp, w[1] - w[3]^2 == 0)
@constraint(nlp, w[2] - 1 - w[4]^2 == 0)
# @constraint(nlp, w[1]^2 + 2*w[2]^2 - 3 - w[4]^2 == 0)
# @constraint(nlp, x[1] >= 0)
# @constraint(nlp, x[1]^2 + 2*x[2]^2 - 3 >= 0)
optimize!(nlp)
solution_summary(nlp)
w_optimal = [value(w[i]) for i ∈ eachindex(w)]

# x_optimal = [value(x[i]) for i ∈ eachindex(x)]
# y_optimal = [value(y[i]) for i ∈ eachindex(y)]
# pr = generate_pr(objective, x0, params=pALP, problemType=problemType; objectiveString=objectiveString)