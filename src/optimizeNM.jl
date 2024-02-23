include("types.jl")
include("nelderMead.jl")

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
    @unpack actions = solverState
    actions[:sort] += 1
    fevals += (n+1)
    Fk = F0
    @pack! solverState = fevals, actions

    solState = SolStateNMType(Xk=X0, Fk=Fk, Delta=Delta0)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k = solverState
        # @show k
        @unpack Xk, Fk = solState

        # saving the current iterates to solState
        Xkm1, Fkm1 = Xk, Fk
        @pack! solState = Xkm1, Fkm1
        
        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_NM = printOrNot & verbose_ls

        Xkp1, fkp1, actions_1NM = nelderMead(Xk, f, pDict, verbose=printOrNot_NM)

        
        @unpack actions = solverState
        # @show actions
        # @show actions_1NM
        actions = merge(+, actions, actions_1NM)
        # @show actions
        @pack! solverState = actions

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = Xkp1[:, 1]
        fvals[k] = fkp1

        Deltak = simplexDiameter(Xk)
        @pack! solState = Deltak

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