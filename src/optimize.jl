# module optimize

# export optimize

function optimize(pr; 
    verbose::Bool=false, 
    log::Bool=true, 
    itrStart::Int64=1)

    # Initial settings
    dftol = pr.alg.dftol
    progress = pr.alg.progress
    maxiter = pr.alg.maxiter
    x0 = pr.x0
    x = x0
    obj = pr.objective
    p = pr.p
    M = max(size(p.data, 1), 1)
    fnext = 1e10
    fₖ = obj(x0, p, getGradientToo=false)
    n = length(x)
    itr = 1
    fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals = zeros(Float64, n, maxiter)
    
    myprintln(true, "Begin with the solver:")
    
    while abs(fnext - fₖ) ≥ dftol && itr ≤ maxiter
        printOrNot = verbose && (itr % progress == 0)
        # printOrNot = false
        myprintln(printOrNot, "Iteration $(itr):", log=true)
        fₖ, ∇fₖ = obj(x, p)
        pₖ = findDirection(pr, ∇fₖ)
        α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, itrStart=itrStart)
        # α, x, fnext, backtrackNum = linesearch_parallel(pr, x, pₖ, itrStart=itrStart)
        fvals[itr] = fnext
        αvals[itr] = α
        backtrackVals[itr] = backtrackNum
        xvals[:, itr] = x
        itr += 1
    end
    
    if itr > maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(itr) iterations 😄"
        myprintln(true, statusMessage)
        # truncating arrays as they weren't filled to capacity
        fvals, αvals, backtrackVals, xvals = [arr[1:itr] for arr in (fvals, αvals, backtrackVals, xvals)]
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
    log::Bool=true)
    
    obj = pr.objective
    p = pr.p
    isStrongWolfe = (pr.alg.linesearch == "StrongWolfe")
    c₁ = pr.alg.c1
    β = 1 / 2^(itrStart-1)
    xnext = copy(xnow)
    fₖ, ∇fₖ = obj(xnow, p)
    fnext = fₖ
    itr_search_for_α = itrStart-1

    while itr_search_for_α ≤ itrMax
        @inbounds xnext .= xnow .+ β .* pₖ
        fnext, ∇fnext = obj(xnext, p)
        comparison_val = fₖ + c₁ * β * dot(∇fₖ, pₖ)

        if fnext ≤ comparison_val
            if isStrongWolfe && abs(dot(∇fnext, pₖ)) < abs(c₁ * dot(∇fₖ, pₖ))
                β /= 2
                itr_search_for_α += 1
            else
                break
            end
        else
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

# function linesearch(pr::NamedTuple, xnow::Vector{Float64}, 
#     pₖ::Vector{Float64};
#     itrMax::Int64=50,
#     itrStart::Int64=1,
#     verbose::Bool=false,
#     log::Bool=true)
    
#     obj = pr.objective
#     p = pr.p
#     linesearch = pr.alg.linesearch
#     c₁ = pr.alg.c1
#     c₂ = pr.alg.c2
#     β = 1/2^(itrStart-1)
#     diff = β*pₖ
#     xnext = xnow+diff
#     fₖ, ∇fₖ = obj(xnow, p, verbose=verbose, log=log)
#     fnext = fₖ
#     itr_search_for_α = itrStart-1
#     myprintln(verbose, "Current value of F, fₖ = $(fₖ)", log=log)
#     armijoSatisfied = false
#     strongWolfeSatisfied = false
#     if linesearch == "StrongWolfe"
#         while !strongWolfeSatisfied && itr_search_for_α ≤ itrMax
#             diff = β*pₖ
#             myprintln(false, "Let's shift x by $(diff)", log=log)
#             xnext = xnow+diff
#             fnext = obj(xnext, p, getGradientToo=false)
#             # println(c₁*β*∇fₖ'*pₖ)
#             myprintln(false, "To be compared against: $(fₖ + c₁*β*∇fₖ'*pₖ)", log=log)
#             if fnext ≤ fₖ + c₁*β*∇fₖ'*pₖ
#                 myprintln(verbose, "Armijo condition satisfied for β = $(β)", log=log)
#                 fnext, ∇fnext = obj(xnext, p)
#                 if abs(∇fnext'*pₖ) ≥ abs(c₂*∇fₖ'*pₖ)
#                     myprintln(verbose, "Curvature condition satisfied for β = $(β)", log=log)
#                     strongWolfeSatisfied = true
#                 else
#                     itr_search_for_α += 1
#                     myprintln(false, "Curvature condition NOT satisfied for β = $(β)", log=log)
#                     β /= 2
#                     myprintln(false, "Line Search Iterations = $(itr_search_for_α)", log=log)
#                 end
#             else
#                 itr_search_for_α += 1
#                 myprintln(false, "Armijo condition NOT satisfied for β = $(β)", log=log)
#                 β /= 2
#                 myprintln(false, "Line Search Iterations = $(itr_search_for_α)", log=log)
#             end 
#         end
#     elseif linesearch == "Armijo"
#         # fₖ, ∇fₖ = obj(xnow, p)
#         while !armijoSatisfied && itr_search_for_α ≤ itrMax
#             diff = β*pₖ
#             myprintln(false, "Let's shift x by $(diff)", log=log)
#             xnext = xnow+diff
#             fnext = obj(xnext, p, getGradientToo=false)
#             # println(c₁*β*∇fₖ'*pₖ)
#             myprintln(false, "To be compared against: $(fₖ + c₁*β*∇fₖ'*pₖ)", log=log)
#             if fnext ≤ fₖ + c₁*β*∇fₖ'*pₖ
#                 myprintln(verbose, "Armijo condition satisfied for β = $(β)", log=log)
#                 armijoSatisfied = true
#             else
#                 itr_search_for_α += 1
#                 myprintln(false, "Armijo condition NOT satisfied for β = $(β)", log=log)
#                 β /= 2
#                 myprintln(false, "Line Search Iterations = $(itr_search_for_α)", log=log)
#             end 
#         end
#     else 
#         @error "Unknown linesearch condition"
#     end
    
#     if itr_search_for_α > itrMax
#         @error "Line Search failed at point x = $(xnext) despite $(itr_search_for_α) iterations."
#     end

#     α = β
#     return (α=α, x=xnext, f=fnext, backtracks=itr_search_for_α) 
# end

# end