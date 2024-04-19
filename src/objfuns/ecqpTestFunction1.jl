include("objective.jl")
include("functionQP.jl")

using Parameters

# Define the quadratic and linear terms of the objective function
G = [2 0; 0 2]
c = [-4; -6]
# Define the constraints
A = [1 1]
b = [1]

# Initial guess
x0 = A\b


pECQP = Dict(:G=>G, :c=>c, :A=>A, :b=>b )

objective = QPObjectiveFunction
objectiveString = "ecqpTestFunction1"
params = pECQP

pr = generate_pr(objective, x0, params=params, problemType="ECQP"; objectiveString=objectiveString)
