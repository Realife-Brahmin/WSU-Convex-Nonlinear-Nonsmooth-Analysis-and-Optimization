include("types.jl")
include("nelderMead.jl")

"""
    optimizeNM(pr; verbose=false, verbose_ls=false, log=true, log_path="./logging/")

Conducts optimization using the Nelder-Mead simplex method, equipped with extensive logging and verbosity for in-depth analysis and debugging. It's designed to tackle complex optimization problems where derivative information is unavailable or unreliable.

# Parameters
- `pr`: A problem definition object that encapsulates the objective function, its parameters, and initial guesses, alongside specific algorithmic configurations like step sizes and termination criteria.
    - `objective`: The objective function to minimize. It should accept a vector of parameters and a dictionary of additional parameters, returning a scalar value.
    - `x0`: Initial guess for the parameters, serving as the starting point for the optimization.
    - `p`: Additional parameters required by the objective function, packed into a dictionary.
    - `alg`: A dictionary containing algorithm-specific parameters such as `alpha`, `beta`, `gamma`, `delta`, `progress`, `maxiter`, and `DeltaTol`.
- `verbose`: If set to `true`, enables the output of detailed information about each iteration's progress and decisions made by the algorithm.
- `verbose_ls`: Controls verbosity for specific aspects of the algorithm, such as line search or simplex adjustments, offering granular insights into the algorithm's execution.
- `log`: Activates the logging of detailed optimization steps and results to a file, aiding in post-optimization analysis.
- `log_path`: Specifies the directory path where the log files will be stored, with each log file named uniquely based on the objective function, algorithm method, and iteration count to prevent overwrite and ensure easy retrieval.

# Operational Details
The function initiates the optimization process by constructing an initial simplex around the provided starting point (`x0`) using a deviation factor to ensure diversity among the simplex vertices. It then enters a loop, iteratively applying Nelder-Mead operations—reflection, expansion, contraction, and shrinkage—guided by the algorithm's parameters (`alpha`, `beta`, `gamma`, `delta`) to explore the parameter space in search of a minimum.

At each iteration, the function evaluates the objective function at simplex vertices, sorts the simplex based on these evaluations, and decides on the next move based on the NM algorithm's logic. It tracks the number of function evaluations, updates the simplex, and checks against termination criteria: exceeding the maximum number of iterations (`maxiter`) or the simplex size falling below a threshold (`DeltaTol`).

# Logging and Debugging
- Detailed logs including the iteration count, current best value, simplex status, and actions taken are generated, which can be directed to the console or to a log file based on `verbose` and `log` flags.
- The logging mechanism is designed to assist in understanding the optimization path, diagnosing potential issues, and fine-tuning algorithm parameters.

# Returns
A comprehensive results dictionary encompassing:
- `converged`: A boolean indicating whether the optimization successfully converged based on the criteria.
- `statusMessage`: A message detailing the outcome of the optimization process.
- `xvals`, `fvals`: Arrays tracking the parameter values and corresponding function evaluations at each iteration.
- `fevals`: The total number of objective function evaluations performed.
- `cause`: Reasons for stopping the optimization, providing insights into the termination condition.
- `pr`, `solState`, `solverState`: Objects encapsulating the problem definition, current solution state, and solver state, respectively, offering a snapshot of the optimization process at termination.

# Examples
To optimize a given function `your_objective_function` starting from an initial guess `[initial_guess]` with specified parameters and algorithm settings:
```julia
pr = Dict(
    :objective => your_objective_function,
    :x0 => [initial_guess],
    :p => Dict(parameter_dictionary),
    :alg => Dict(:method => "Nelder-Mead", :maxiter => 100, :DeltaTol => 1e-5)
)
result = optimizeNM(pr, verbose=true, log=true)
```
"""
function optimizeNM(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg[:method]*"_"*string(pr.alg[:maxiter])*".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    
    solverState = SolverStateNMType()

    @unpack alpha, beta, gamma, delta = pr.alg

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    DeltaTol = pr.alg[:DeltaTol]

    x0 = pr.x0
    n = length(x0)
    xk = x0

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, n, maxiter)

    # doing this even though NM requires multiple f evals
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    f = pr.objective
    pDict = pr.p

    fk = f(x0, pDict, getGradientToo=false)
    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    X00 = createInitialSimplexFromOnePoint(x0, deviationFactor=0.1) # this simplex is currently unsorted

    X0, F0 = sortSimplex(X00, f, pDict)
    Delta0 = simplexDiameter(X0)
    @unpack actions, fevals = solverState
    actions[:sort] += 1
    fevals += (n+1)
    f0 = F0[:, 1]
    @pack! solverState = actions, fevals

    solState = SolStateNMType(Xk=X0, fk=f0, Deltak=Delta0)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k = solverState
        @unpack Xk, fk = solState

        # saving the current iterates to solState
        Xkm1, fkm1 = Xk, fk
        @pack! solState = Xkm1, fkm1

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_NM = printOrNot & verbose_ls

        Xkp1, fkp1, actions_1NM, fevals_1NM = nelderMead(Xk, f, pDict, verbose=printOrNot_NM)

        @unpack actions, fevals = solverState

        fevals += fevals_1NM
        actions = merge(+, actions, actions_1NM)
        @pack! solverState = actions, fevals

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = Xkp1[:, 1]
        fvals[k] = fkp1

        Deltak = simplexDiameter(Xk)
        @pack! solState = Deltak # pre-emptively packing it into the solState, as it won't be mutated

        Xk, fk = Xkp1, fkp1
        @pack! solState = Xk, fk

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
        elseif Deltak < DeltaTol
            push!(causeForStopping, "Simplex size lower limit reached!")
            keepIterationsGoing = false
        end

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

    res = (converged=converged, statusMessage=statusMessage,    xvals=xvals, fvals=fvals, fevals=fevals, cause=causeForStopping, pr=pr, solState=solState,
    solverState=solverState)

    res = trim_array(res, k - 1)

    return res
    
end