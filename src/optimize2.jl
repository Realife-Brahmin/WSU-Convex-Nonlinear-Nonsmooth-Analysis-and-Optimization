include("linesearches.jl");
include("findDirection.jl");
include("types.jl");

include("getCandidateStep.jl");
include("checkQualityOfCandidateStep.jl");
include("updateTRegionParams.jl");
include("updateTRModelParams.jl");

function optimize2(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    method = pr.alg[:method]
    linesearchMethod = pr.alg[:linesearch]
    if method == "TrustRegion"
        linesearchMethod = "QuasiNewton-SR1"
    end

    log_txt = log_path*"log_"*string(pr.objective)*"_"*method*"_"*linesearchMethod*"_"*string(pr.alg[:maxiter])*".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    solverState = SolverStateType()
    solState = SolStateType(xk=pr.x0)

    # Initial settings
    fevals = 0
    gevals = 0
    dftol = pr.alg[:dftol]
    gtol = pr.alg[:gtol]
    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    x0 = pr.x0
    xk = x0

    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    obj = pr.objective
    p = pr.p
    data = p[:data]
    M = max(size(data, 1), 1)

    fk = obj(x0, p, getGradientToo=false)
    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)
    @pack! solState = fk

    if method == "QuasiNewton"
        QNState = QNStateType()
    elseif method == "ConjugateGradientDescent"
        CGState = CGStateType()
    elseif method == "TrustRegion"
        SR1params = SR1paramsType()
        TRparams = TRparamsType()
    end
    
    fevals += 1
    @pack! solverState = fevals

    n = length(xk)

    fvals, αvals, gmagvals = [zeros(Float64, maxiter) for _ in 1:3]
    backtrackVals = zeros(Int64, maxiter)
    xvals, gvals = [zeros(Float64, n, maxiter) for _ in 1:2]
    
    myprintln(true, "Begin with the solver:", log=log, log_path=log_txt)
    keepIterationsGoing = true
    causeForStopping = []

    justRestarted = false # automatically false if not doing CGD, and if doing CGD and latest β was not zero.

    while keepIterationsGoing

        @unpack k = solverState

        printOrNot = verbose && ( (k - 1) % progress == 0)
        printOrNot_ls = printOrNot & verbose_ls


        myprintln(printOrNot, "Iteration $(k):", log_path=log_txt)

        fk, gk = obj(xk, p)
        @checkForNaN fk
        @checkForNaN gk

        gmagk = sum(abs.(gk))
        usingCGD = false
        fevals += 1
        gevals += 1

        @pack! solState = fk, gk, gmagk
        @pack! solverState = fevals, gevals

        if method == "QuasiNewton"
            @pack! QNState = k, xk, fk, gk
            pk, QNState = findDirection(pr, gk, QNState=QNState)
            # @show pk

        elseif method == "ConjugateGradientDescent"
            usingCGD = true
            @pack! CGState = k, gk
            pk, CGState = findDirection(pr, gk, CGState=CGState)
            @unpack justRestarted = CGState 
        
        elseif method == "TrustRegion"
            SR1params = updateTRModelParams(SR1params, solState)

        else
            pk = findDirection(pr, gk)
        end

        if method != "TrustRegion"
        
            @pack! solState = pk 

            if linesearchMethod == "StrongWolfe"

                solState, solverState = StrongWolfe(pr, solState, solverState,
                verbose=printOrNot_ls)


            elseif linesearchMethod == "Armijo"
                @error "Armijo no longer supported."
            
            else
                @error "Unknown linesearch method"
            end

            @unpack success_ls = solverState
            if ~success_ls
                myprintln(true, "Line search failed... Bad direction or optimal point?")
                push!(causeForStopping, "LineSearch failed.")
                keepIterationsGoing = false
            end

        elseif method == "TrustRegion"

            @unpack xk = solState
            y = xk

            keepFindingCandidate = true
            
            while keepFindingCandidate
                pj = getCandidateStep(SR1params, TRparams)
                ρk, fj, solverState = checkQualityOfCandidateStep(SR1params, pr, solState, pj, solverState)
                
                y, TRparams = updateTRegionParams(TRparams, ρk, pj, y)

                @unpack Delta, Delta_min = TRparams
                @unpack xk = solState
                
                if y != xk
                    keepFindingCandidate = false
                    myprintln(verbose_ls, "A new valid candidate step has been found")
                    xkm1 = xk
                    xk = y
                    @pack! solState = xkm1, xk

                elseif Delta < Delta_min
                    keepFindingCandidate = false
                    myprintln(verbose_ls, "Trust Region Radius too small now.")
                else
                    myprintln(verbose_ls, "Keep Finding Candidate")
                end

            end

        else
            @error "floc"
        end

        @unpack xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk = solState
        # @show xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk
        myprintln(printOrNot, "Iteration $(k): x = $(xk) is a better point with new fval = $(fk).", log_path=log_txt)

        if !usingCGD && !justRestarted && abs(fk - fkm1) < dftol
            push!(causeForStopping, "Barely changing fval")
            keepIterationsGoing = false
        end
        if !usingCGD && !justRestarted && gmagkm1 < gtol
            push!(causeForStopping, "Too small gradient at previous step.")
            keepIterationsGoing = false
        end
        if !justRestarted && gmagk < gtol
            push!(causeForStopping, "Too small gradient at latest step.")
            keepIterationsGoing = false
        end
        if method == "TrustRegion"
            @unpack Delta, Delta_min = TRparams
            if Delta < Delta_min
                push!(causeForStopping, "Trust Region Radius too small")
                keepIterationsGoing = false
            end
        end
        if k == maxiter
            push!(causeForStopping, "Too many iterations")
            keepIterationsGoing = false
        end

        @unpack Hk, alphak = solState
        @unpack alpha_evals = solverState

        fvals[k] = fk
        αvals[k] = alphak
        gvals[:, k] = gk
        gmagvals[k] = gmagk
        backtrackVals[k] = alpha_evals
        xvals[:, k] = xk

        k += 1

        @pack! solverState = k
        @pack! solState = k

    end
    
    @unpack k = solverState

    if k ≥ maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log=log,  log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(k) iterations 😄"
        myprintln(true, statusMessage, log=log, log_path=log_txt)
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, 
    αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, gmagvals=gmagvals, 
    gvals=gvals, M=M, fevals=fevals, gevals=gevals, cause=causeForStopping, 
    pr=pr)

    res = trim_array(res, k-1)
    return res
end