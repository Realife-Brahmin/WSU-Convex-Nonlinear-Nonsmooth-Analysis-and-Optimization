
include("augmentedLagrangian.jl")

function optimizeALP(pr;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    objString = pr.objectiveString
    log_txt = log_path * "log_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    solverState = SolverStateALPType()

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    etol = pr.alg[:etol]
    itol = pr.alg[:itol]

    x0 = pr.x0
    n = length(x0)
    xk = x0

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, n, maxiter)

    # doing this even though GA requires multiple f evals
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)

    f = pr.objective
    pALP = pr.p
    @unpack G, c, mE, Ae, be, mI, A, b = pALP
    

    f0 = f(x0, pALP, getGradientToo=false)
    fk = f0
    solState = SolStateALPType(x0, Ae, fk=f0, Wk0=Wk0, itol=itol)
    # @show solState[:Wk]
    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    keepIterationsGoing = true
    causeForStopping = []

    myprintln(verbose, "*"^25, log_path=log_txt)

    while keepIterationsGoing

        @unpack k = solverState

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ALP = printOrNot & verbose_ls

        myprintln(verbose, "-"^10, log_path=log_txt)

        myprintln(printOrNot_ALP, "Iteration k = $(k)", log_path=log_txt)

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
            break
        end

        @unpack xk, fk, Wk = solState
        myprintln(printOrNot_ALP, "Current Working Set: $(Wk)")
        # Since this iteration will be proceeded with, saving the current iterates to solState as the 'previous' iteration values
        km1, xkm1, fkm1, Wkm1 = k, xk, fk, Wk
        @pack! solState = km1, xkm1, fkm1, Wkm1

        myprintln(printOrNot_ALP, "Let's try taking a step saisfying current active constraints.")
        pk = getECQPStep(pr, solState, verbose=printOrNot_ALP, verbose_ls=printOrNot_ALP)

        if norm(pk) < itol # stationary point wrt Wk

            myprintln(printOrNot_ALP, "Step is too small, checking if this is a stationary point wrt all constraints.")

            lambdas = computeLagrangianMultipliersQP(xk, G, c, Awk)

            lambda_min, jI_min = findmin(lambdas)

        elseif norm(pk) >= itol

            myprintln(printOrNot_ALP, "Step pk of size $(norm(pk)) exceeds tolerance, so proceeding with checking if any 'outside' constraints are blocking it.")
        else

            @error("floc")

        end

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = xkp1
        fvals[k] = fkp1 # also incorrect

        xk, fk, Wk = xkp1, fkp1, Wkp1

        @pack! solState = xk, fk, Wk
        @pack! solState = k
        @pack! solverState = k

    end

    @unpack k = solverState

    myprintln(true, "*"^50, log=log, log_path=log_txt)

    if k ≥ maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log=log, log_path=log_txt, color=:red)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(k) iterations 😄"
        myprintln(true, statusMessage, log=log, log_path=log_txt, color=:green)
    end

    @unpack fevals = solverState

    res = (converged=converged, statusMessage=statusMessage, xvals=xvals, fvals=fvals, fevals=fevals, cause=causeForStopping, pr=pr, solState=solState,
        solverState=solverState)

    res = trim_array(res, k)

    return res

end