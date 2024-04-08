include("objective.jl")
include("equalityConstrainedQP.jl")

using Parameters

# Nocedal and Wright Example 1 in QP formulation
# Define the quadratic and linear terms of the objective function
G = [6 2 1;2 5 2; 1 2 4]
c = [-8, -3, -3]
# Define the constraints
A = [1 0 1; 0 1 1]
b = [3, 0]

# Initial guess
x0 = A\b

pECQP = Dict(:G=>G, :c=>c, :A=>A, :b=>b )

objective = equalityConstrainedQP
objectiveString = "ecqpTestFunction2"
params = pECQP

pr = generate_pr(objective, x0, params=params, problemType="ECQP"; objectiveString=objectiveString)
