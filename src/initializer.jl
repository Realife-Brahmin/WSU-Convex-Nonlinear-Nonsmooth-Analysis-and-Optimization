module initializer

using DataFrames

export estimate_omega

function estimate_omega(df::DataFrame)::Float64
    # Extract the time and value columns
    times = df.t
    values = df.V

    # Initialize variables
    peak_times = Float64[]
    increasing = false

    for i in 2:length(values)
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

end