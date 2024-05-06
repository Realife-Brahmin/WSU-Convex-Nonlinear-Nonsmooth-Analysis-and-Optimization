include("objective.jl")
include("../AugmentedLagrangian.jl")

using Parameters

x0 = rand(2)
# x0 = ones(2)
xk = x0
nx = length(x0)

function alpTestFunction01(x, p;
    getGradientToo::Bool=true)

    if length(x) != 2
        @error "alpTestFunction01 expects a length 2 vector"
    end

    f = (x[1] - 1)^2 + (x[2] - 2)^2
    if !getGradientToo
        return f
    elseif getGradientToo
        g = [2 * (x[1] - 1), 2 * (x[2] - 2)]
        return f, g
    else
        @error("floc")
    end

    @error("floc")
end

function cE01(x, p;
    getGradientToo::Bool=true)

    if length(x) != 2
        @error "cE01 expects a length 2 vector"
    end

    mE = 1
    cE = zeros(mE)
    cE[1] = x[1]^2 + x[2]^2 - 1

    if !getGradientToo
        return cE
    elseif getGradientToo
        nx = length(x)
        hE = zeros(mE, nx)
        hE[1, 1] = 2*x[1]
        hE[1, 2] = 2*x[2]
        return cE, hE
    else
        @error("floc")
    end

    @error("floc")

end

mE = 1

function cI01(x, p;
    getGradientToo::Bool=true)

    if length(x) != 2
        @error "cI01 expects a length 2 vector"
    end

    mI = 2
    cI = zeros(mI)
    cI[1] = x[1]^3 + x[2]
    cI[2] = x[1]^2 + 2*x[2]^2 - 1.1
    if !getGradientToo
        return cI
    elseif getGradientToo
        nx = length(x)
        hI = zeros(mI, nx)
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
y0 = rand(mI)
# y0 = zeros(mI)
yk = y0
m = mE+mI
w0 = vcat(x0, y0)
wk = w0

objective = alpTestFunction01
objectiveOriginal = alpTestFunction01
objectiveString = "alpTestFunction01"
problemType = "Constrained"
econ = cE01
icon = cI01
pDictALP = Dict(:nx=>nx, :m=>m, :mE=>mE, :econ=>econ, :mI=>mI, :icon=>icon)
pr = generate_pr(objective, w0, params=pDictALP, problemType=problemType; objectiveString=objectiveString)

# objectiveUnc = ALOBJ
# subroutineCall = true
# muk = 1.0
# m = mE+mI
# lambdak = ones(m)
# addendum = Dict(:subroutineCall => subroutineCall, :lambda => lambdak, :mu => muk, :objective => objective, :objectiveUnc => objectiveUnc)
# pDictUnc = merge(deepcopy(pDictALP), addendum)
# # xk = x0
# F, G = ALOBJ(wk, pDictUnc, getGradientToo=true)

# using JuMP, Ipopt

# model = Model(Ipopt.Optimizer)
# # x0 = [-1.2, 1.0]
# # x0 = [1, 2]
# slackifyInequalities = false
# slackifyInequalities = true
# # nx = length(x0)
# @variable(model, x[i=1:nx], start = x0[i])
# @variable(model, y[i=1:mI], start = y0[i])
# @objective(model, Min, (x[1] - 1)^2 + (x[2] - 2)^2)
# @constraint(model, x[1]^2 + x[2]^2 - 1 == 0) # regular equality constraint

# if slackifyInequalities
#     @constraint(model, x[1]^3 + x[2] - y[1]^2 == 0) 
#     @constraint(model, x[1]^2 + 2 * x[2]^2 - 1.1 - y[2]^2 == 0)
# else
#     @constraint(model, x[1]^3 + x[2] >= 0)
#     @constraint(model, x[1]^2 + 2 * x[2]^2 - 1.1 >= 0)
# end

# optimize!(model)
# xopt = [value(x[i]) for i ∈ 1:nx]
# yopt = [value(y[i]) for i ∈ 1:mI]
# f_optimal = objective_value(model)
# println("Optimal Variables x: ", xopt)
# println("Optimal Slack Variables y: ", yopt)
# println("Optimal objective value: ", f_optimal)

