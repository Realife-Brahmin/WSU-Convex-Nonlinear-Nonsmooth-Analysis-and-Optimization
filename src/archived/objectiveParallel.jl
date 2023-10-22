function dampedSHM_Parallel(x::Vector{Float64}, 
    # p::NamedTuple{(:df, :params), Tuple{DataFrame, Vector{Float64}}};
    p::NamedTuple{(:data, :params), Tuple{Matrix{Float64}, Vector{Float64}}};
    getGradientToo::Bool=true,
    verbose::Bool=false,
    log::Bool=true)

    # df = p.df
    # y = df.y
    # t = df.x
    # M = length(y)
    # df = p.df
    data = p.data
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf
    n = length(x)
    A₀, A, τ, ω, α, ϕ = x
    f = 0.0;
    g = zeros(Float64, n)
    A₀, A, τ, ω, α, ϕ = x
    
    f_atomic = Threads.Atomic{Float64}(0.0)  # Make f atomic
    g_atomic = [Threads.Atomic{Float64}(0.0) for _ in 1:n]  # Make each component of g atomic

    @threads for k = 1:M
        tₖ = t[k]
        yₖ = y[k]
        expₖ = exp(-tₖ/τ)
        Sₖ = expₖ*sin((ω+α*tₖ)tₖ + ϕ)
        Cₖ = expₖ*cos((ω+α*tₖ)tₖ + ϕ)

        ŷₖ = A₀ + A*Sₖ
        Δyₖ = ŷₖ - yₖ
        
        Threads.atomic_add!(f_atomic, (1/M)*Δyₖ^2)  # Use atomic add for f
        if getGradientToo
            Δg = (2/M)* Δyₖ * [1, Sₖ, A*tₖ*(τ^-2)*Sₖ, A*tₖ*Cₖ, A*tₖ^2*Cₖ, A*Cₖ]
            for j = 1:n
                Threads.atomic_add!(g_atomic[j], Δg[j])  # Use atomic add for each component of g
            end
        end
    end

    f = f_atomic[]  # Extract the value from atomic
    g = [g_atomic[j][] for j = 1:n]  # Extract values from atomic array

    if getGradientToo
        return f, g
    else
        return f
    end
end