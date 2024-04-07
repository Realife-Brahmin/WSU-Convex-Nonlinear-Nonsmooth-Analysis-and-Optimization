include("objective.jl")
include("equalityConstrainedQP.jl")

using Parameters

function ecqpTestFunction1(x::Vector{Float64},
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    @unpack G, c = pDict[:params]
    f = 1 // 2 * transpose(x) * G * x + transpose(x) * c

    if !getGradientToo
        return f
    elseif getGradientToo
        g = G*x + c
        return f, g
    else
        @error "floc"
    end

end

# Define the quadratic and linear terms of the objective function
G = [2 0; 0 2]
c = [-4; -6]
# c = zeros(2)
# Define the constraints
A = [1 1]
b = [1]

# Initial guess
x0 = A\b
# x0 = [0.5; 0.5]
# frac = rand()
# x0 = [frac, 1-frac]
# x0 = [1.0, 0.0]

pECQP = Dict(:G=>G, :c=>c, :A=>A, :b=>b )

objective = equalityConstrainedQP
objectiveString = "ecqpTestFunction1"
params = pECQP

pr = generate_pr(objective, x0, params=params, problemType="ECQP"; objectiveString=objectiveString)

# f0 = equalityConstrainedQP(w0, pECQP)

# solState = SolStatePGCGType(w0, G, c, A, fk=f0)
