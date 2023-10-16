function linesearchArmijo(pr::NamedTuple, xk::Vector{Float64}, 
    pₖ::Vector{Float64};
    ρ = 0.5,
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")
    
    fevals_ls = 0
    gevals_ls = 0
    obj = pr.objective
    p = pr.p
    c₁ = pr.alg.c1
    β = ρ^(itrStart-1)
    
    xnext = copy(xk)
    if pr.alg.method == "GradientDescent"
        fₖ = obj(xk, p, getGradientToo=false)
        ∇fₖ = -pₖ
    else
        fₖ, ∇fₖ = obj(xk, p)
        gevals_ls += 1
    end
    fevals_ls += 1
    
    fnext = fₖ
    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    itr_search_for_α = itrStart-1

    while itr_search_for_α ≤ itrMax
        xnext = xk + β * pₖ
        myprintln(verbose, "Let's try shifting x to $(xnext)", log_path=log_txt)
        
        comparison_val = fₖ + c₁ * β * pₖ'*∇fₖ
        fnext = obj(xnext, p, getGradientToo=false)
        fevals_ls += 1
        
        if fnext ≤ comparison_val
            myprintln(verbose, "Armijo condition satisfied for β = $(β)", log_path=log_txt)
            break;
        else
            myprintln(false, "Armijo condition NOT satisfied for β = $(β)", log=log, log_path=log_txt)
            β *= ρ
            itr_search_for_α += 1
        end
    end

    if itr_search_for_α > itrMax
        @error "Line Search failed at point x = $(xnext) despite $(itr_search_for_α) iterations."
    end

    α = β
    return (α=α, x=xnext, f=fnext, backtracks=itr_search_for_α, fevals=fevals_ls, gevals=gevals_ls) 
end

function strongWolfeBisection(pr::NamedTuple, xk, pk; 
    α_min=0, α_max=100, 
    itrMax=50, tol=1e-8,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    c₁ = pr.alg.c1
    c₂ = pr.alg.c2
    obj = pr.objective
    p = pr.p
    ϕ0, ϕprime0 = ((f, g) -> (f, pk'*g))(obj(xk, p)...)
    fevals_ls = 1
    gevals_ls = 1
    # ϕprime0 = pk'*g(xk)
    iteration = 0
    xmid = xk
    ϕ_mid = ϕ0
    α_mid = 0
    strongWolfeSatisfied = false
    while iteration < itrMax && abs(α_max - α_min) > tol
        α_mid = 0.5 * (α_min + α_max)
        xmid = xk + α_mid * pk
        myprintln(verbose, "Let's try shifting x to $(xmid)", log_path=log_txt)
        # Compute function value and gradient at α_mid
        ϕ_mid = obj(xmid, p, getGradientToo=false)
        fevals_ls += 1
        # Check the Armijo condition
        if ϕ_mid > ϕ0 + c₁ * α_mid * ϕprime0
            α_max = α_mid
        else
            myprintln(verbose, "StrongWolfe1 condition satisfied for α = $(α_mid)", log_path=log_txt)
            # ϕprime_mid = pk'*g(xmid)
            ϕmid, ϕprime_mid = ((f, g) -> (f, pk'*g))(obj(xmid, p)...)
            fevals_ls += 1
            gevals_ls += 1
            # Check the curvature condition
            if abs(ϕprime_mid) <= c₂ * abs(ϕprime0)
                myprintln(verbose, "StrongWolfe2 condition satisfied for α = $(α_mid)", log_path=log_txt)
                strongWolfeSatisfied = true
                return (α=α_mid, x=xmid, f=ϕ_mid, backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls)
            else
                myprintln(verbose, "StrongWolfe2 condition NOT satisfied for α = $(α_mid)", log_path=log_txt)
                α_min = α_mid
            end
        end
        
        iteration += 1
    end
    
    # if abs(α_max - α_min) < tol
    #     @warn "Strong Wolfe conditions not met even after $iteration iterations, but tolerance of gap of $tol between αmax and αmin achieved." 
    # else
    #     @warn "Strong Wolfe line search did not converge in $itrMax iterations, Returning best α with tolerance gap of $(abs(α_max - α_min)) between αmax and αmin"
    # end
    return (α=α_mid, x=xmid, f=ϕ_mid, backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls)
end

# α_min = 0.0
# α_max = 100.0

# ϕ(x) = sum((x - π*ones(Float64, length(x))).^4)
# dϕ(x) = 4*(x-π*ones(Float64, length(x))).^3
# xk = [0.0]
# xkp1, fkp1, αk = strongWolfeBisection(ϕ, dϕ, xk, -dϕ(xk), α_min, α_max, ϕ(xk),  -dϕ(xk)'*dϕ(xk))
# println("x_k+1 = $xkp1,  f_k+1 = $fkp1,  α_k = $αk")
# # plot(f, α_min, α_max, label="f(x)", xlabel="x", ylabel="f(x)", title="Plot of f(x)")

# f(x) = sum((x - ones(Float64, length(x))).^2)
# g(x) = 2 * (x - ones(Float64, length(x)))
# xk = [2.0, 2.0]
# xkp1, fkp1, αk,  = strongWolfeBisection(f, g, xk, -g(xk), α_min, α_max, f(xk), -g(xk)'*g(xk))
# println("x_k+1 = $xkp1,  f_k+1 = $fkp1,  α_k = $αk")



# f(x) = sum((x-4*ones(Float64, length(x))).^2)
# g(x) = 2(x-4*ones(Float64, length(x)))
# xk = [0]
# xkp1, fkp1, αk = strongWolfeBisection(f, g, xk, -g(xk), α_min, α_max, f(xk),  -g(xk)'*g(xk))
# xkp1, fkp1, αk = strongWolfeBisection(f, g, xkp1, -g(xkp1), α_min, α_max, f(xkp1),  -g(xkp1)'*g(xkp1))

# println("x_k+1 = $xkp1,  f_k+1 = $fkp1,  α_k = $αk")
# # plot(f, α_min, α_max, label="f(x)", xlabel="x", ylabel="f(x)", title="Plot of f(x)")

# f(x) = sum(x.^4 - 8*x.^3 + 18*x.^2)
# g(x) = 4*x.^3 - 24*x.^2 + 36*x
# xk = [5]
# xkp1, fkp1, αk = strongWolfeBisection(f, g, xk, -g(xk), α_min, α_max, f(xk),  -g(xk)'*g(xk))
# println("x_k+1 = $xkp1,  f_k+1 = $fkp1,  α_k = $αk")
# xkp1, fkp1, αk = strongWolfeBisection(f, g, xkp1, -g(xkp1), α_min, α_max, f(xkp1),  -g(xkp1)'*g(xkp1))

