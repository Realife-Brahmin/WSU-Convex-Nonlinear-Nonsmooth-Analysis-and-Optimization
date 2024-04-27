include("objective.jl")

using Parameters

# Define the quadratic and linear terms of the objective function
G = [2 0; 0 2]
c = [-4; -6]
# Define the constraints
Ae = [1 1]
be = [1]

# Initial guess
x0 = Ae \ be


pECQP = Dict(:G => G, :c => c, :Ae => Ae, :be => be)

objective = QPObjectiveFunction
objectiveString = "ecqpTestFunction1"
params = pECQP

pr = generate_pr(objective, x0, params=params, problemType="ECQP"; objectiveString=objectiveString)