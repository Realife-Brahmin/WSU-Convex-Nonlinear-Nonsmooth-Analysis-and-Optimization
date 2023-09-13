# module initializer

using DataFrames

# export estimate_omega

function estimate_x0(df::DataFrame, 
    x::Vector{Num})::Vector{Float64}

    x0 = zeros(Float64, length(x))
    ω = estimate_omega(df)
    α = 0
    A₀ = estimate_A0(df)
    A = estimate_A(df)
    x0[1] = A₀
    x0[2] = A
    x0[3] = 1
    x0[4] = ω
    x0[5] = α # α decay is almost zero
    x0[6] = estimate_phi(df, A₀, A,  ω)

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

# function estimate_phi(df, ω)
#     # Find the first zero crossing after some transient time (for simplicity, let's say after the first 5% of the data)
#     start_idx = round(Int, 0.05 * length(df.t))
#     zero_crossing_time = -1

#     # Search for the zero crossing by looking where the sign changes
#     for i in start_idx:(length(df.t) - 1)
#         if df.V[i] * df.V[i+1] < 0  # A change in sign indicates a zero crossing
#             zero_crossing_time = df.t[i] + (df.t[i+1] - df.t[i]) * abs(df.V[i]) / (abs(df.V[i]) + abs(df.V[i+1]))
#             break
#         end
#     end

#     # Expected time of zero crossing for unshifted sine wave would be when ωt = kπ
#     expected_time = π / ω
    
#     # The difference in time gives the phase offset
#     φ = ω * (expected_time - zero_crossing_time)
    
#     # Make sure φ is in the interval [0, 2π]
#     φ = mod(φ, 2π)
    
#     return φ
# end

function estimate_phi(df::DataFrame, A₀::Float64, A::Float64, ω::Float64, t_peak::Float64)
    adjusted_V = df.V .- A₀
    V_peak = A * exp(-t_peak / τ)  # Based on the model without the sine term
    # Calculate the desired value at t_peak using the adjusted V
    desired_value_at_peak = adjusted_V[df.t .≈ t_peak][1] / V_peak

    # Make sure the value is within the range for arcsine
    if abs(desired_value_at_peak) > 1.0
        println("Warning: Desired value at peak is out of domain for arcsine. Adjusting to lie within [-1, 1].")
        desired_value_at_peak = sign(desired_value_at_peak)
    end

    φ = asin(desired_value_at_peak) - ω * t_peak
    return φ
end

# end