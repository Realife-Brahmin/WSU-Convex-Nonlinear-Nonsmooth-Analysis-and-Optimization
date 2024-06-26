include("geneticAlgorithm.jl")
# include("nelderMead.jl")
include("sampleSpace.jl")
include("types.jl")

"""
    optimizeGA(pr; verbose=false, verbose_ls=false, log=true, log_path="./logging/")

Executes a genetic algorithm (GA) optimization process based on predefined parameters and configurations. It iteratively evolves a population towards better solutions, assessing convergence criteria and logging progress throughout the optimization journey.

## Parameters:
- `pr`: A structure containing the problem definition, including the objective function, algorithm parameters, and initial conditions.
- `verbose`: Enables general logging of the optimization process when set to `true`.
- `verbose_ls`: Enables detailed logging specifically for the genetic algorithm's operations.
- `log`: If `true`, enables logging to a file specified by `log_path`.
- `log_path`: Directory path where the log file is saved, structured based on the objective function name, method, and iteration limit.

## Returns:
- A tuple containing the final state of the optimization process, including convergence status, final population and fitness values, total number of function evaluations, and a record of actions taken during the optimization.

## Methodology:
1. Initializes logging and solver state, unpacks genetic algorithm parameters from `pr`.
2. Generates an initial population and evaluates its fitness.
3. Enters the main optimization loop, performing parent selection, crossover, mutation, and population update in each iteration until convergence criteria are met.
4. Logs detailed progress and decisions at each step if `verbose` options are enabled.
5. Assesses the optimization's success based on fitness improvement and iteration count, returning comprehensive results and diagnostics.

## Examples:
```julia
# Define the problem structure `pr` with the objective function, algorithm parameters, and initial conditions
pr = defineProblem(objective=myObjectiveFunction, alg=setupGAParameters(), x0=initialCondition)

# Execute the GA optimization
results = optimizeGA(pr, verbose=true, log=true)
```
"""
function optimizeGA(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg[:method]*"_"*string(pr.alg[:maxiter])*".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    
    solverState = SolverStateGAType()

    # @unpack genetic parameters = pr.alg

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    fvalRepeatTol = pr.alg[:fvalRepeatTol]
    popSize = pr.alg[:popSize]
    deviation = pr.alg[:deviation]
    delta = pr.alg[:delta]
    Dist = pr.alg[:Dist]
    parentsSurvive = pr.alg[:parentsSurvive]
    dftol = pr.alg[:dftol]

    x0 = pr.x0
    n = length(x0)
    xk = x0

    popSize = min(popSize, n+1)

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, n, maxiter)

    # doing this even though GA requires multiple f evals
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    f = pr.objective
    pDict = pr.p

    fk = f(x0, pDict, getGradientToo=false)
    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    X0, F0, f0 = createInitialPopulation(x0, popSize, f, pDict, deviationFactor=deviation, verbose=verbose) # the first element of the next generation is expected to be the best one

    @unpack fevals = solverState
    fevals += popSize
    @pack! solverState = fevals

    solState = SolStateGAType(Xk=X0, Fk=F0, fk=f0)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k, fvalRepeats = solverState

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
            break
        elseif fvalRepeats >= fvalRepeatTol
            push!(causeForStopping, "Fitness no longer generationally improving!")
            keepIterationsGoing = false
            break
        end

        @unpack Xk, Fk, fk = solState
        # saving the current iterates to solState
        Xkm1, Fkm1, fkm1 = Xk, Fk, fk
        @pack! solState = Xkm1, Fkm1, fkm1

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_GA = printOrNot & verbose_ls

        myprintln(printOrNot, "Iteration k = $(k)", log_path=log_txt)

        Xkp1, Fkp1, fkp1, fevals_1GA, actions_1GA = deriveNextGeneration(Xk, Fk, fk, f, pDict, delta=delta, deviation=deviation, Dist=Dist,
        parentsSurvive = parentsSurvive, dftol = dftol, verbose=printOrNot_GA) # first element should be the best one

        @unpack actions, fevals = solverState

        fevals += fevals_1GA
        actions = merge(+, actions, actions_1GA)
        @pack! solverState = actions, fevals

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = Xkp1[:, 1]
        fvals[k] = fkp1 # also incorrect
        
        @unpack actions = solverState
        if abs(fkp1 - fk) < dftol
            fvalRepeats += 1
            actions[:genFitnessNotImproved] += 1
            myprintln(printOrNot_GA, "Generation Fitness not improved.", log_path=log_txt)
        else
            fvalRepeats = 0
            actions[:genFitnessImproved] += 1
            myprintln(printOrNot_GA, "Fittest individual now even fitter!", log_path=log_txt)
        end
        @pack! solverState = actions, fvalRepeats

        @pack! solState = fvalRepeats # pre-emptively packing it into the solState, as it won't be mutated

        Xk, Fk, fk = Xkp1, Fkp1, fkp1
        @pack! solState = Xk, Fk, fk


        @pack! solState = k
        @pack! solverState = k

    end

    @unpack k = solverState

    if k ≥ maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log=log, log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(k) iterations 😄"
        myprintln(true, statusMessage, log=log, log_path=log_txt)
    end

    @unpack fevals = solverState

    res = (converged=converged, statusMessage=statusMessage, xvals=xvals, fvals=fvals, fevals=fevals, cause=causeForStopping, pr=pr, solState=solState,
    solverState=solverState)

    res = trim_array(res, k - 1)

    return res
    
end