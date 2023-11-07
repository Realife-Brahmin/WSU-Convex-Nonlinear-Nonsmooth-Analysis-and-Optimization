function StrongWolfeBisection(pr::NamedTuple, xk, pk; 
    α_min=0, α_max=1, 
    itrMax=50, tol=1e-10,
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
    iteration = 0
    xmid = xk
    ϕ_mid = ϕ0
    ϕprime_mid = 0
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

            ϕmid, ϕprime_mid = ((f, g) -> (f, pk'*g))(obj(xmid, p)...)
            fevals_ls += 1
            gevals_ls += 1
            # Check the curvature condition
            if abs(ϕprime_mid) <= c₂ * abs(ϕprime0)
                myprintln(verbose, "StrongWolfe2 condition satisfied for α = $(α_mid)", log_path=log_txt)
                strongWolfeSatisfied = true

                ans = (α=α_mid, x=xmid, f=ϕ_mid, gmag=abs(ϕprime_mid), backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls, success=true)

                return ans
            else

                myprintln(verbose, "StrongWolfe2 condition NOT satisfied for α = $(α_mid)", log_path=log_txt)
                α_min = α_mid
            end
        end
        
        iteration += 1
    end
    
    if abs(α_max - α_min) ≤ tol && iteration ≤ itrMax
        @warn "Strong Wolfe line search did not converge in $(iteration) iterations, Returning best α with tolerance gap of $(abs(α_max - α_min)) between αmax and αmin"
        success=false
    elseif abs(α_max - α_min) > tol && iteration == itrMax
        @warn "Strong Wolfe line search did not converge in $(itrMax) iterations, Returning best α with tolerance gap of $(abs(α_max - α_min)) between αmax and αmin"
        success=false
    else
        @error "Bad condition"
        println(α_max - α_min)
        println(iteration)
    end
    
    ans = (α=α_mid, x=xmid, f=ϕ_mid, gmag=abs(ϕprime_mid), backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls, success=false)

    return ans

    @error "floc"

end


