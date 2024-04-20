include("objective.jl")
include("functionQP.jl")
include("../activeSetQP.jl")

using Parameters

G = [2 0;0 2]
c = [-2, -5]
c = Float64.(c)
n = length(c)
c0 = sum([1, 2.5].^2) # constant offset to the QP objective
lb = myfill(c, -Inf)
ub = myfill(c, Inf)

Ae = zeros(0, n)
be = zeros(0)
A = [1 -2;-1 -2;-1 2;1 0;0 1]
b = [-2, -6, -2, 0, 0]
mE, mI = length(be), length(b)

Wk0 = [3, 5]
pQP = @packDict "{G, c, lb, ub, c0, mE, Ae, be, mI, A, b, Wk0}"

# w0 = computeFeasiblePointForLinearConstraints(pQP)
w0 = [2, 0]
objective = QPObjectiveFunction
objectiveOriginal = QPObjectiveFunction
objectiveString = "asqpTestFunction1" # Ex 16.3 in Nocedal and Wright
# params = pQP

pr = generate_pr(objective, w0, params=pQP, problemType="QP"; objectiveString=objectiveString)