include("objective.jl")

using Parameters

# function cE01(x, p, getGradientToo::Bool=false)
    

#     if !getGradientToo
#         return f
#     elseif getGradientToo
#         return f, g
#     else
#         @error("floc")
#     end

#     @error("floc")

# end

# function cI01(x, p, getGradientToo::Bool=false)


#     if !getGradientToo
#         return f
#     elseif getGradientToo
#         return f, g
#     else
#         @error("floc")
#     end

#     @error("floc")

# end

G = [2 0; 0 2]
c = [-2, -5]
c = Float64.(c)
n = length(c)
c0 = sum([1, 2.5] .^ 2) # constant offset to the QP objective


Ae = zeros(0, n)
be = zeros(0)
A = [1 -2; -1 -2; -1 2; 1 0; 0 1]
b = [-2, -6, -2, 0, 0]
mE, mI = length(be), length(b)

# pQP = @packDict "{G, c, lb, ub, c0, mE, Ae, be, mI, A, b, Wk0}"


# objective = QPObjectiveFunction
# objectiveOriginal = QPObjectiveFunction
objectiveString = "alpTestFunction01" # Ex 16.3 in Nocedal and Wright
# params = pQP

# pr = generate_pr(objective, w0, params=pQP, problemType="ALP"; objectiveString=objectiveString)

using JuMP, NLPModelsJuMP, Percival
nlp = Model(NLPModelsJuMP.Optimizer)
set_attribute(nlp, "solver", Percival.PercivalSolver)
@variable(nlp, x[i=1:2], start = [-1.2; 1.0][i])
@objective(nlp, Min, (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2)
@constraint(nlp, x[1]^2 + x[2]^2 == 1)
optimize!(nlp)
solution_summary(nlp)
# value.(nlp)
x_optimal = [value(x[i]) for i ∈ eachindex(x)]s