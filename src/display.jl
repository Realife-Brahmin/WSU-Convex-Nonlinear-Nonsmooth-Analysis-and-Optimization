
using DataFrames
using LaTeXStrings
using Parameters
using Plots
using Printf

# include("utilities.jl")
include("objective.jl")

# Extract fields from the 'res' structure and assign them to local variables.
# The purpose of this macro is to avoid repetitive assignments.
# Usage: @unpack_vars res status fvals αvals backtrackVals xvals
macro unpack_vars(struct_var, fields...)
    quote
        # For each specified field, assign the field's value from the structure to a local variable.
        $(Expr(:block, [:(local $(field) = $(struct_var).$(field)) for field in fields]...))
    end
end

function scatter_voltage_vs_time(df::DataFrame)
    gr()  # Ensure you're using the GR backend for Plots.jl which supports LaTeX rendering

    theme(:dark)  # Setting a dark theme
    
    # Convert time to milliseconds
    time_ms = df.t .* 1000

    # Create the scatter plot
    scatter(time_ms, df.V, 
        label=L"Voltage $(V)$",  # Using LaTeXStrings with \text for regular text
        xlabel=L"Time $(ms)$",  # Using LaTeXStrings with \text for regular text
        ylabel=L"Voltage $(mV)$",  # Using LaTeXStrings with \text for regular text
        title="Voltage vs Time",
        markersize=4,
        color=:purple,
        alpha=0.8
    )
    
    # Additional customization if needed
    # background_color_legend=:transparent
end

function showresults(res::NamedTuple)
    @unpack converged, statusMessage, fvals, xvals, backtrackVals = res
    nnztol = 1e-8 # variable having a value below this will NOT be counted in the list of non-zero variables. For printing purposes only. 
    n = size(xvals, 1)
    itr = size(xvals, 2)
    x☆ = xvals[:, itr]
    minBacktracks, maxBacktracks = extrema(backtrackVals)

    v = true
    myprintln(v, "****************************")
    myprintln(v, statusMessage)
    if converged == true
        fvalPrefix = "Optimal MSE value = "
        xvalMessage = "Nonzero Optimal Variables:"
    elseif converged == false
        fvalPrefix = "Best MSE value = "
        xvalMessage = "Nonzero Best Known Variables:"
    else
        @error "bad condition"
    end
    fvalMessage = fvalPrefix*"$(fvals[itr])"
    myprintln(v, fvalMessage)
    myprintln(v, "***************************")
    myprintln(v, xvalMessage)
    for k in 1:n
        if abs(x☆[k]) > nnztol     
            myprintln(v, "x☆[$(k)] = $(x☆[k])")
        end
    end
    myprintln(v, "***************************")
end


using LaTeXStrings
using Plots

function plotresults(pr::NamedTuple, res::NamedTuple; savePlot::Bool=true, saveLocation::String="processedData/")
    @unpack converged, statusMessage, fvals, xvals, backtrackVals = res
    n = size(xvals, 1)
    itr = size(xvals, 2)
    x☆ = xvals[:, itr]

    # Convert time to milliseconds
    time_ms = pr.df.t .* 1000
    
    # Calculate dampedSHM values using optimal parameters
    predicted_values = [pr.objective(x☆, t, getGradientToo=false) for t in pr.df.t]

    # Plot
    gr()  # Ensure you're using the GR backend for Plots.jl which supports LaTeX rendering
    theme(:dark)  # Setting a dark theme

    p = scatter(time_ms, pr.df.V, 
        label=L"Original Data $(V)$",  # Using LaTeXStrings for original data label
        xlabel=L"Time $(ms)$",  # Using LaTeXStrings 
        ylabel=L"Voltage $(mV)$",  # Using LaTeXStrings 
        title="Voltage vs Time",
        markersize=4,
        color=:purple,
        alpha=0.8
    )
    
    # Plot dampedSHM predictions
    plot!(p, time_ms, predicted_values, 
        label="Predicted with dampedSHM", 
        linewidth=2,
        color=:cyan,
        alpha=0.8
    )

    # Save the plot if savePlot is true
    if savePlot
        # Create the filename using the provided logic
        filename = saveLocation * string(pr.objective) * "_" * pr.alg.method * "_" * pr.alg.linesearch * "_$(@sprintf("%.0e", itr)).png"
        
        # Ensure directory exists
        mkpath(dirname(filename))
        
        # Save the plot to the generated filename
        savefig(p, filename)
    end

    # Display plot
    display(p)
end

