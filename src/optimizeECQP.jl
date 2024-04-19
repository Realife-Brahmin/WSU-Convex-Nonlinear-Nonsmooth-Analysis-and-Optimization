include("projectedGradientConjugateGradient.jl")

"""
    optimizeECQP(pr; verbose=false, verbose_ls=false, log=true, log_path="./logging/") -> Dict

Solves an Extended Constrained Quadratic Programming (ECQP) problem using the Projected Gradient Conjugate Gradient (PGCG) method, providing detailed control over logging and verbosity levels throughout the optimization process.

# Arguments
- `pr`: A problem representation `NamedTuple` or `Dict` containing all necessary information for the ECQP problem, including the objective function, constraints, and algorithm-specific parameters.

# Keyword Arguments
- `verbose::Bool=false`: Enables detailed printing of the optimization process if set to `true`.
- `verbose_ls::Bool=false`: Controls the verbosity of line search operations within the PGCG method. Effective when `verbose` is also `true`.
- `log::Bool=true`: Activates logging of the optimization process to a file.
- `log_path::String="./logging/"`: Specifies the directory path for saving log files related to the optimization process.

# Returns
- `Dict`: A dictionary containing the optimization results and states, including convergence status, the final values of variables, function evaluations, and solver states.

# Notes
- The function initializes solver and solution states, managing iterations through the PGCG method to minimize the objective function subject to equality constraints.
- Logging is handled conditionally based on `log` and `verbose` flags, enabling both console output and file-based logging of the optimization journey.
- The solver dynamically adapts to convergence criteria and iteration limits, ensuring robust handling of the ECQP problem.
- `pr.alg` must be properly defined in the input `pr`, including keys for algorithmic parameters like `method`, `maxiter`, `dftol`, and `etol`.

# Example
```julia
# Define ECQP problem parameters and algorithmic settings
G = [2 0; 0 2]
A = [1 1]
c = [-1; -1]
b = [0]
x0 = [0.5; 0.5]
alg_settings = Dict(:method => "PGCG", :maxiter => 100, :dftol => 1e-5, :etol => 1e-6, :progress => 10)

pr = (objectiveString="MyObjective", alg=alg_settings, x0=x0, p=Dict(:params => Dict(:G => G, :c => c, :A => A, :b => b)))

# Execute the optimization
results = optimizeECQP(pr, verbose=true, log=true)
```
"""
function optimizeECQP(pr;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    objString = pr.objectiveString
    log_txt = log_path * "log_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run


    solverState = SolverStateECQPType()

    # @unpack genetic parameters = pr.alg

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    dftol = pr.alg[:dftol]
    etol = pr.alg[:etol]

    x0 = pr.x0
    n = length(x0)
    xk = x0

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, n, maxiter)

    # doing this even though GA requires multiple f evals
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    f = pr.objective
    pECQP = pr.p
    @unpack G, c, A, b = pECQP

    AAT = A*transpose(A)

    f0 = f(x0, pECQP, getGradientToo=false)
    fk = f0
    solState = SolStatePGCGType(x0, G, c, A, fk=f0, etol=etol)

    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k = solverState

        @unpack xk, gk, dk, rk = solState
        # saving the current iterates to solState
        km1, xkm1, gkm1, dkm1, rkm1 = k, xk, gk, dk, rk
        num = transpose(rk)*gk

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
            break
        elseif num < etol
            push!(causeForStopping, "Convergence Tolerance Reached")
            keepIterationsGoing = false
            break
        end

        @pack! solState = km1, xkm1, gkm1, dkm1, rkm1

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ECQP = printOrNot & verbose_ls

        myprintln(printOrNot, "Iteration k = $(k)", log_path=log_txt)

        @unpack actions, fevals = solverState

        xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG = solveForNextPGCGIterate(xk, gk, dk, rk, num, G, A, AAT, verbose=printOrNot_ECQP) # last two oargs are bogus

        fkp1 = f(xkp1, pECQP, getGradientToo=false)
        fevals += 1
        fevals += fevals_1PGCG # bogus, does nothing
        actions = merge(+, actions, actions_1PGCG) # bogus, does nothing

        @pack! solverState = actions, fevals

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = xkp1
        fvals[k] = fkp1 # also incorrect

        xk, gk, dk, rk = xkp1, gkp1, dkp1, rkp1

        @pack! solState = xk, gk, dk, rk

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

    res = trim_array(res, k)

    return res

end