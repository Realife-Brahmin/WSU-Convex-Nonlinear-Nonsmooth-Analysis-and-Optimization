include("helperFunctions.jl")

function deriveNextGeneration(Xk::Matrix, f::Function, pDict::Dict;
    verbose::Bool = false)

    n, p = size(Xk)
    Xkp1 = myfill(Xk, -75.0)
    Fkp1 = myzeros(Fk)
    survived = zeros(p)
    popSize = 0
    keepAddingToPopulation = true
    fevals = 0
    actions = Dict(
        :bothParentsSurvived => 0, :onlyOneParentSurvived => 0, :noParentSurvived => 0, :crossover => 0, :genFitnessImproved => 0,
        :genFitnessNotImproved => 0, :mutation => 0,
        :mutationFailure => 0, :mutationSuccess => 0, :OneChild => 0, :OneChildOneParent => 0,
        :OneChildBothParents => 0, :parentsSelected => 0)


    # add the top survivor, it is assumed that Xk, Fk have the best values at the front
    Xkp1[:, 1] = Xk[:, 1]
    Fkp1[1] = Fk[1]
    popSize += 1
    # TODO increment suitable actions for the top survivor
    
    while keepAddingToPopulation

        p1, p2 = selectParents(Xk, Fk)
        actions[:parentsSelected] += 1

        o = crossover(p1, p2) # offspring
        actions[:crossover] += 1

        om = mutate(o, f, pDict)
        

    end

    return Xkp1, Fkp1, fevals_1GA, actions_1GA
end