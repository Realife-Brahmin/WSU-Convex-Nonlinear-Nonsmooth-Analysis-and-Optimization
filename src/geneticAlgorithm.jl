include("helperFunctions.jl")

function deriveNextGeneration(Xk::Matrix, f::Function, pDict::Dict;
    verbose::Bool=false)

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
    popSize += 1
    Xkp1[:, popSize] = Xk[:, 1]
    Fkp1[popSize] = Fk[1]

    # TODO increment suitable actions for the top survivor

    while keepAddingToPopulation

        if popSize >= p
            keepAddingToPopulation = false
            myprintln(verbose, "Next Generation Complete.")
        end

        p1, p2 = selectParents(Xk, Fk)
        actions[:parentsSelected] += 1

        o = crossover(p1, p2) # offspring
        actions[:crossover] += 1

        om, mutations_1mut = mutate(o, f, pDict)
        fevals += 1
        actions[:mutations] += mutations_1mut

        Xkp1, Fkp1, popSize, survived, actions_1dAA = decideAndAdd(p1, p2, om, Xkp1, Fkp1, popSize, pDict) # select the child, and if choosing to let the parent survive by default, the best min(2, p-popSize) parents
        actions = merge(+, actions, actions_1dAA)

    end

    Xkp1, Fkp1 = fittestFirst(Xkp1, Fkp1)

    if Fkp1[1] > Fk[1]
        actions[:genFitnessImproved] += 1
        myprintln(verbose, "Fittest individual now even fitter.")
    else
        actions[:genFitnessNotImproved] += 1
        myprintln(verbose, "Fittest individual has same fitness as previous generation.")
    end

    return Xkp1, Fkp1, fevals, actions
end