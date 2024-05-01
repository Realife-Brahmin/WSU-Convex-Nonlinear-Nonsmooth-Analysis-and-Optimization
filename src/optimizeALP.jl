
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
    dxtol = pr.alg[:dxtol]
    gtol = pr.alg[:gtol]
    mutol = pr.alg[:mutol]

    w0 = pr.x0 # actually w0
    N = length(w0)
    wk = w0

    f = pr.objective
    pALP = pr.p
    @unpack n, m, mE, econ, mI, icon = pALP

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, N, maxiter)

    myprintln(verbose, "Starting with initial point w = $(wk).", log_path=log_txt)



    lambdak0 = rand(m)
    x0, y0 = w0[1:n], w0[n+1:end]
    f0 = f(x0, pALP, getGradientToo=false)
    fk = f0
    solState = SolStateALPType(w0, lambdak=lambdak0, fk=f0, etol=etol) # again, x0 is actually w0

    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    keepIterationsGoing = true
    causeForStopping = []

    myprintln(verbose, "*"^25, log_path=log_txt)

    while keepIterationsGoing

        @unpack k = solverState
        muk = solState[:muk] # separately calling it here to not break my pattern later
        
        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ALP = printOrNot & verbose_ls

        myprintln(verbose, "-"^10, log_path=log_txt)

        myprintln(printOrNot_ALP, "Iteration k = $(k)", log_path=log_txt)

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
            break

        elseif muk >= mutol

            myprintln(printOrNot_ALP, "μ has exceeded the limit of $(mutol) !")
            keepIterationsGoing = false
            break

        end

        @unpack wk, fk, gk, lambdak, muk, tk = solState

        # Since this iteration will be proceeded with, saving the current iterates to solState as the 'previous' iteration values
        km1, wkm1, fkm1, gkm1, lambdakm1, mukm1, tkm1 = k, wk, fk, gk, lambdak, muk, tk
        @pack! solState = km1, wkm1, fkm1, gkm1, lambdakm1, mukm1, tkm1

        myprintln(printOrNot_ALP, "Let's try solving for the Augmented Lagrangian problem with λ = $(lambdak) and μ = $(muk).")
        wkp1, fkp1, cxkp1, iter_unc = solveAugmentedLagrangianFunction(pr, solState, verbose=printOrNot_ALP, verbose_ls=printOrNot_ALP)
        
        normdxk = norm(wkp1 - wk)

        if normdxk < dxtol

            myprintln(printOrNot_ALP, "Change is too small, perhaps I will check if equality constraints are satisfied.")
            lambdakp1 = lambdak
            mukp1 = muk
            tkp1 = tk
            push!(causeForStopping, "Change in decision variables too small!")
            keepIterationsGoing = false
            break

        elseif normdxk >= dxtol

            myprintln(printOrNot_ALP, "Change of size $(normdxk) exceeds tolerance, so will update λ and μ.")
            lambdakp1 = lambdak - transpose(muk)*cxkp1
            mukp1 = muk*( 1 + 10*exp(-iter_unc/n) )
            tkp1 = tk/2

        else

            @error("floc")

        end

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = wkp1
        fvals[k] = fkp1

        wk, fk, lambdak, muk, tk = wkp1, fkp1, lambdakp1, mukp1, tkp1

        @pack! solState = wk, fk, gk, lambdak, muk, tk
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

    xopt, fopt = extractBestResults(pr, k, xvals, fvals)

    @unpack fevals = solverState

    res = (converged=converged, statusMessage=statusMessage, 
            iter=k,
            xvals=xvals, xopt=xopt, 
            fvals=fvals, fopt=fopt, 
            fevals=fevals,
            cause=causeForStopping, 
            pr=pr, 
            solState=solState,
            solverState=solverState)

    res = trim_array(res, k)

    return res

end