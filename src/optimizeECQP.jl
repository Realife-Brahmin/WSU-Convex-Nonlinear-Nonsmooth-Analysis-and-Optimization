include("projectedGradientConjugateGradient.jl")

"""
    optimizeECQP(pr; verbose=false, verbose_ls=false, log=true, log_path="./logging/")

Executes an optimization cycle using the Extended Conjugate Quadratic Programming (ECQP) method.

# Arguments
- `pr`: The problem record containing all necessary settings, including the initial point, problem data, and algorithm configuration.

# Keyword Arguments
- `verbose::Bool=false`: Enables verbose output, providing detailed logs of the optimization process.
- `verbose_ls::Bool=false`: Enables verbose output specifically for line search operations.
- `log::Bool=true`: Flag to enable or disable logging of the optimization process.
- `log_path::String="./logging/"`: The directory path where the log files will be stored.

# Returns
- A tuple containing:
    - `converged::Bool`: A boolean flag indicating if the optimization converged.
    - `statusMessage::String`: A message describing the outcome of the optimization.
    - `xvals::Matrix{Float64}`: The history of x values throughout the iterations.
    - `fvals::Vector{Float64}`: The history of objective function values throughout the iterations.
    - `fevals::Int`: Total number of function evaluations performed.
    - Other state and control data from the optimization process.

# Description
This function orchestrates the ECQP solver by initializing the solver state, managing iterations, and handling convergence checks. It leverages active set strategies to efficiently navigate the solution space and applies the Projected Gradient Conjugate Gradient (PGCG) method within each iteration to adjust the search direction and step size. 

Verbose outputs and logging are managed based on user settings, providing insights into the internal workings and progression of the solver. The function is designed to handle both standalone and subroutine calls within a larger optimization framework, such as Active Set Quadratic Programming (ASQP).

# Notes
- The function performs a pre-run check to ensure the feasibility of the initial point against the equality constraints `Ae * x = be`.
- If discrepancies are found that exceed the error tolerance (`etol`), the solver will not proceed, ensuring robustness in the face of potentially erroneous inputs.

# Examples
```julia
# Define the problem settings and initial conditions
pr = {
    objectiveString: "QuadraticObjective",
    p: {G: Matrix, c: Vector, Ae: Matrix, be: Vector},
    alg: {method: "PGCG", maxiter: 100, progress: 10, etol: 1e-8},
    x0: [initial_guess]
}

# Perform optimization
result = optimizeECQP(pr, verbose=true, log_path="/path/to/logs/")

# Analyze the results
println("Optimization converged: ", result.converged)
println("Final objective value: ", result.fvals[end])
```
"""
function optimizeECQP(pr;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    objString = pr.objectiveString
    pECQP = pr.p

    if !haskey(pECQP, :subroutineCall)
        myprintln(false, "Calling ECQP Solver Independently.")
        subroutineCall = false
        @unpack G, c, Ae, be = pECQP
    else
        myprintln(false, "Calling ECQP Solver as a subroutine for ASQP.")
        @unpack G, c, Ae, be, subroutineCall = pECQP
    end

    verbose = verbose && !subroutineCall

    log_txt = log_path * "log_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"
    if !subroutineCall
        if isfile(log_txt)
            rm(log_txt)
        end # remove logfile if present for the run
    end

    solverState = SolverStateECQPType()

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    etol = pr.alg[:etol]

    x0 = pr.x0
    n = length(x0)
    xk = x0
    
    myprintln(verbose, "Starting with initial point x = $(xk).", log=verbose, log_path=log_txt)

    myprintln(false, "But first, a quick safety check to make sure we're inserting a feasible xk into our ECQP: What's Ae*xk - be?")
    discrepancy = Ae*xk - be
    if norm(discrepancy) < etol
        myprintln(false, "It's a zero vector as it should be. All good!")
    else
        error("The residual is non-zero!. The ECQP solver cannot be allowed to proceed!")
    end

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, n, maxiter)

    f = pr.objective


    AAT = Ae*transpose(Ae)

    f0 = f(x0, pECQP, getGradientToo=false)
    fk = f0
    solState = SolStatePGCGType(x0, G, c, Ae, fk=f0, etol=etol)

    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log=verbose, log_path=log_txt)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k = solverState

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ECQP = printOrNot & verbose_ls

        myprintln(printOrNot, "Iteration k = $(k)", log=verbose, log_path=log_txt)

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

        @unpack actions, fevals = solverState

        xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG = solveForNextPGCGIterate(xk, gk, dk, rk, num, G, Ae, AAT, verbose=printOrNot_ECQP) # last two oargs are bogus

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
        myprintln(!subroutineCall, statusMessage, log=log, log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(k) iterations 😄"
        myprintln(!subroutineCall, statusMessage, log=log, log_path=log_txt)
    end

    @unpack fevals = solverState

    res = (converged=converged, statusMessage=statusMessage, xvals=xvals, fvals=fvals, fevals=fevals, cause=causeForStopping, pr=pr, solState=solState,
        solverState=solverState)

    res = trim_array(res, k)

    return res

end