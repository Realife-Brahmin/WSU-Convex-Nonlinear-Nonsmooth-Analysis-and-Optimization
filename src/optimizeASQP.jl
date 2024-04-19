
include("activeSetQP.jl")

function optimizeASQP(pr;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    objString = pr.objectiveString
    log_txt = log_path * "log_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    solverState = SolverStateASQPType()

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
    pASQP = pr.p
    @unpack G, c, mE, Ae, be, mI, A, b = pASQP

    f0 = f(x0, pASQP, getGradientToo=false)
    fk = f0
    solState = SolStateASQPType(x0, Ae, fk=f0, itol=itol)

    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k = solverState
        myprintln(printOrNot, "Iteration k = $(k)", log_path=log_txt)

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ASQP = printOrNot & verbose_ls

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
            break
        end

        @unpack xk, fk, Wk = solState
        # Since this iteration will be proceeded with, saving the current iterates to solState as the 'previous' iteration values
        km1, xkm1, fkm1, Wkm1 = k, xk, fk, Wk
        @pack! solState = km1, xkm1, fkm1, Wkm1

        pk = getStepDirection() # write a function to invoke ECQP
        Iall = collect(mE+1:mE+mI) # vector of indices for all inequality constraints [4, 5, 6, 7, 8] where mE = 3 mI = 5
        WIk = Wk[mE+1:end] # contains only indices for inequality constraints (like [5, 7, 8])
        notWIk = setdiff(Iall, WIk) # [4, 6]
        Awk = A[WIk .- mE, :] # contains only working set inequality constraints (like A[ [2, 4, 5], :])
        Anotwk = A[notWIk .- mE, :] # A[ [1, ]]
        bwk = A[WIk .- mE] # A[ [2, 4, 5], :]
        bnotwk = b[notWIk .- mE]
        if norm(pk) < tol # stationary point wrt Wk

            lambdas = computeLagrangianMultipliers(xk, G, c, Awk)

            lambda_min, jI_min = findmin(lambdas)

            if lambda_min >= 0
                myprintln(verbose, "0 vector is the only feasible descent direction")

                push!(causeForStopping, "No feasible improvement step possible from here.")
                keepIterationsGoing = false
                break

            else
                j_min = jI_min + mE
                myprintln(verbose, "Working set constraint $(j_min) prevents a feasible step the hardest.")
                myprintln(verbose, "So removing constraint $(j_min) from the working set")
                xkp1 = xk
                fkp1 = fk
                Wkp1 = vcat(Wk[1:j_min-1], Wk[j_min+1:end])

            end

        elseif norm(pk) > tol

            alphak = 1.0
            jBlockClosest = 0

            for idx ∈ eachindex(Anotwk)
                den = Anotwk[idx, :]*pk
                if den < 0
                    distance = (bnotwk[idx] - Anotwk[idx, :]*xk)/den
                    if distance < alphak
                        jBlockClosest = notWIk[idx]
                        myprintln(printOrNot_ASQP, "Constraint $(jBlockClosest) outside the working set is blocking the current step the earliest among the checked ones.")
                        alphak = distance
                    end
                end
            end

            if jBlockClosest == 0
                myprintln(printOrNot_ASQP, "No inequality constraint outside the working set is blocking our current step.")
                Wkp1 = Wk
                xkp1 = xk + pk
                fkp1 = f(xkp1, pASQP, getGradientToo=false)

            elseif jBlockClosest > 0
                myprintln(printOrNot_ASQP, "Constraint $(jBlockClosest) is blocking the step first, so we'll restrict our step accordingly and add it to the working set.")
                Wkp1 = sort!(push!(Wk, jBlockClosest))
                xkp1 = xk + alphak*pk
                fkp1 = f(xkp1, pASQP, getGradientToo=false)

            else
                @error("floc")
            
            end

        else

            @error("floc")

        end


        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = xkp1
        fvals[k] = fkp1 # also incorrect

        @pack! solState = xk, fk, Wk
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