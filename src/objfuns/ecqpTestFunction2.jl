include("objective.jl")
include("functionQP.jl")

using LinearAlgebra
using Parameters

# Nocedal and Wright Example 1 in QP formulation
# Define the quadratic and linear terms of the objective function
G = [6 2 1;2 5 2; 1 2 4]
c = [-8, -3, -3]
# Define the constraints
Ae = [1 0 1; 0 1 1]
be = [3, 0]

# Initial guess
p = vec(nullspace(Ae))

x00 = Ae\be
x0 = x00 + rand()*p
pECQP = Dict(:G=>G, :c=>c, :Ae=>Ae, :be=>be )

objective = QPObjectiveFunction
objectiveString = "ecqpTestFunction2"
params = pECQP

pr = generate_pr(objective, x0, params=params, problemType="ECQP"; objectiveString=objectiveString)
