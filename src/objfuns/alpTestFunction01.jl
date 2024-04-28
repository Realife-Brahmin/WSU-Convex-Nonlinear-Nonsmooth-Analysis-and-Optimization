include("objective.jl")

using Parameters

function cE01(x, p, getGradientToo::Bool=false)
    
    f = x[1]^2 + x[2]^2 - 1
    if !getGradientToo
        return f
    elseif getGradientToo
        g = [2*x[1], 2*x[2]]
        return f, g
    else
        @error("floc")
    end

    @error("floc")

end

function cI01(x, p, getGradientToo::Bool=false)
    # x[1]^3 + x[2] >= 0
    f = x[1]^3 + x[2]
    if !getGradientToo
        return f
    elseif getGradientToo
        g = [3*x[1]^2, 1]
        return f, g
    else
        @error("floc")
    end

    @error("floc")

end

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
pALP = Dict(:mE=>1, :econ=>cE01, :mI=>1, :icon=>cI01)
params = pALP

x0 = [-1.2, 1.0]
using JuMP, NLPModelsJuMP, Percival
nlp = Model(NLPModelsJuMP.Optimizer)
set_attribute(nlp, "solver", Percival.PercivalSolver)
@variable(nlp, x[i=1:2], start = x0[i])
@objective(nlp, Min, (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2)
@constraint(nlp, x[1]^2 + x[2]^2 == 1)
@constraint(nlp, x[1]^3 + x[2] >= 0)
optimize!(nlp)
solution_summary(nlp)
# value.(nlp)
x_optimal = [value(x[i]) for i ∈ eachindex(x)]

# pr = generate_pr(objective, x0, params=pALP, problemType=problemType; objectiveString=objectiveString)