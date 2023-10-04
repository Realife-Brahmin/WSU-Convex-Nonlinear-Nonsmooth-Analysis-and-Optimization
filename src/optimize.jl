# module optimize

# export optimize

function optimize(pr; 
    verbose::Bool=false, 
    log::Bool=true,
    log_path::String="./logging/",
    itrStart::Int64=1)

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run
    # Initial settings
    dftol = pr.alg.dftol
    progress = pr.alg.progress
    maxiter = pr.alg.maxiter
    x0 = pr.x0
    x = x0

    myprintln(verbose, "Starting with initial point x = $(x).", log_path=log_txt)
    obj = pr.objective
    @show p = pr.p
    M = max(size(p.data, 1), 1)
    fnext = 1e10
    fₖ = obj(x0, p, getGradientToo=false)
    n = length(x)
    itr = 1
    fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals = zeros(Float64, n, maxiter)
    
    myprintln(true, "Begin with the solver:", log_path=log_txt)
    
    while abs(fnext - fₖ) ≥ dftol && itr ≤ maxiter
        printOrNot = verbose && (itr % progress == 0)
        # printOrNot = false
        myprintln(printOrNot, "Iteration $(itr):", log_path=log_txt)
        fₖ, ∇fₖ = obj(x, p)
        pₖ = findDirection(pr, ∇fₖ)
        α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, itrStart=itrStart, verbose=printOrNot)
        myprintln(printOrNot, "Iteration $(itr): x = $(x) is a better point with new fval = $(fnext).", log_path=log_txt)
        fvals[itr] = fnext
        αvals[itr] = α
        backtrackVals[itr] = backtrackNum
        xvals[:, itr] = x
        itr += 1
    end
    
    if itr > maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(itr) iterations 😄"
        myprintln(true, statusMessage, log_path=log_txt)
        # truncating arrays as they weren't filled to capacity
        fvals, αvals, backtrackVals, xvals = [arr[1:itr-1] for arr in (fvals, αvals, backtrackVals, xvals)]
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, M=M)

    return res
end


function findDirection(pr::NamedTuple, ∇fnow::Vector{Float64};
    verbose::Bool=false)::Vector{Float64}
    method = pr.alg.method
    n = length(∇fnow)
    if method == "GradientDescent"
        # Bₖ = I(n)
        # pₖ = -Bₖ*∇fnow
        pₖ = -∇fnow
    elseif method == "ConjugateGradientDescent"
        @error "Currently not formulated for this method"
    elseif method == "QuasiNewton"
        @error "Currently not formulated for this method"
    else
        @error "Currently not formulated for this method"
    end

    return pₖ
end

function linesearch(pr::NamedTuple, xnow::Vector{Float64}, 
    pₖ::Vector{Float64};
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")
    
    obj = pr.objective
    p = pr.p
    isStrongWolfe = (pr.alg.linesearch == "StrongWolfe")
    c₁ = pr.alg.c1
    β = 1 / 2^(itrStart-1)
    xnext = copy(xnow)
    fₖ, ∇fₖ = obj(xnow, p)
    fnext = fₖ
    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    itr_search_for_α = itrStart-1

    while itr_search_for_α ≤ itrMax
        xnext .= xnow .+ β .* pₖ
        myprintln(verbose, "Let's try shifting x to $(xnext)", log_path=log_txt)
        fnext, ∇fnext = obj(xnext, p)
        comparison_val = fₖ + c₁ * β * dot(∇fₖ, pₖ)

        if fnext ≤ comparison_val
            myprintln(verbose, "Armijo condition satisfied for β = $(β)", log_path=log_txt)
            if isStrongWolfe && abs(dot(∇fnext, pₖ)) < abs(c₁ * dot(∇fₖ, pₖ))
                myprintln(false, "Curvature condition NOT satisfied for β = $(β)", log_path=log_txt)
                β /= 2
                itr_search_for_α += 1
            else
                break
            end
        else
            myprintln(false, "Armijo condition NOT satisfied for β = $(β)", log=log)
            β /= 2
            itr_search_for_α += 1
        end
    end

    if itr_search_for_α > itrMax
        @error "Line Search failed at point x = $(xnext) despite $(itr_search_for_α) iterations."
    end

    α = β
    return (α=α, x=xnext, f=fnext, backtracks=itr_search_for_α) 
end

# end