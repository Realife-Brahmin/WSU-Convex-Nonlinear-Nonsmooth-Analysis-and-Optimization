# module initializer

using DataFrames
using Symbolics

# export estimate_omega

function estimate_x0(df::DataFrame, 
    x::Vector{Num})::Vector{Float64}

    x0 = zeros(Float64, length(x))
    ω = estimate_omega(df)
    α = 0
    A₀ = estimate_A0(df)
    A = estimate_A(df)
    ϕ = estimate_phi(df, A₀, A,  ω)
    τ = estimate_tau(df, A₀)
    x0[1] = A₀
    x0[2] = A
    x0[3] = τ
    x0[4] = ω
    x0[5] = α # α decay is almost zero
    x0[6] = ϕ

    return x0
end

function estimate_omega(df::DataFrame)::Float64
    # Extract the time and value columns
    times = df.t
    values = df.V

    # Initialize variables
    peak_times = Float64[]
    increasing = false

    for i in 2:lastindex(values)
        if increasing
            if values[i] < values[i-1]
                # A peak is detected
                push!(peak_times, times[i-1])
                increasing = false
            end
        else
            if values[i] > values[i-1]
                increasing = true
            end
        end
    end

    # We need at least two peaks to estimate the frequency
    if length(peak_times) < 2
        error("Couldn't detect sufficient peaks to estimate omega.")
    end

    # Calculate the average period between peaks
    T = mean(diff(peak_times))

    # Calculate omega
    ω = 2 * π / T
    return ω
end

function estimate_A0(df::DataFrame)::Float64
    A₀₀ = mean(df.V)
    return A₀₀
end

function estimate_A(df::DataFrame)::Float64
    # Assume the baseline A₀ is the mean of the voltage data
    A₀ = mean(df.V)
    
    # Adjust the data by subtracting the baseline
    adjusted_V = df.V .- A₀
    
    # Estimate A as the maximum absolute value of the adjusted data
    A = maximum(abs.(adjusted_V))
    
    return A
end

function estimate_phi(df::DataFrame, A₀::Float64, A::Float64, ω::Float64)
    # adjusted_V = df.V .- A₀

    # # Compute the first zero-crossing point after an upward slope
    # # We find where the adjusted voltage crosses zero and has a positive slope
    # zero_crossings = findall(x -> x > 0, diff(sign.(adjusted_V)))
    
    # if isempty(zero_crossings)
    #     throw(ErrorException("Unable to determine zero-crossing from data."))
    # end

    # t_zero_crossing = df.t[zero_crossings[1]]

    # # The phase shift φ would then be determined by comparing this t_zero_crossing to a standard sine function.
    # # Specifically, for a standard sine function, sin(ω*t + φ) = 0 implies that ω*t + φ = 0 (mod π).
    # # Hence, φ = -ω*t_zero_crossing (mod π).

    # φ = mod(-ω * t_zero_crossing, π)
    ϕ = asin((df.V[1] - A₀)/A) - df.t[1]*ω
    return ϕ
end

using Statistics
using DataFrames

function estimate_tau(df::DataFrame, A0_est::Float64)
    # Adjust V values based on the estimated A0
    V_adj = df.V .- A0_est
    
    # Only consider positive adjusted V values for log calculations
    mask = V_adj .> 0
    V_adj_log = log.(V_adj[mask])
    t_values = df.t[mask]
    
    # Compute slope of ln(V_adj) vs t
    slope = cov(t_values, V_adj_log) / var(t_values)
    
    # Return estimated tau
    return -1/slope
end



# end