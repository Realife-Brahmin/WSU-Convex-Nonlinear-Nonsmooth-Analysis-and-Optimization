
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
    time_ms = df.x .* 1000

    # Create the scatter plot
    scatter(time_ms, df.y, 
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

    formatted_iterations = @sprintf("%.0e", itr)

    title_content = "Voltage vs Time\nMethod: $(pr.alg.method)\nLineSearch: $(pr.alg.linesearch)\nIterations: $formatted_iterations\nError (MSE): $(fvals[itr]) mV"

    p = scatter(time_ms, pr.df.V, 
        label=L"Original Data $(V)$",  # Using LaTeXStrings for original data label
        xlabel=L"Time $(ms)$",  # Using LaTeXStrings 
        ylabel=L"Voltage $(mV)$",  # Using LaTeXStrings 
        title="Voltage vs Time",
        # title=title_content,
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

# Define formatting variables
font_size = 8
font_color = :white
font_family = :courier

# Calculate positions for annotations based on the data
max_x = maximum(pr.df.t)
max_y = minimum(pr.df.V)
annotation_spacing = (maximum(pr.df.V) - max_y) * 0.05

# Calculate starting positions for the annotations
start_y = max_y + 4 * annotation_spacing

# Add annotations for method details with monospace font
annotate!(p, [(max_x * 0.95, start_y, text("Method: $(pr.alg.method)", font_color, :right, font_size, font_family)),
             (max_x * 0.95, start_y - annotation_spacing, text("LineSearch: $(pr.alg.linesearch)", font_color, :right, font_size, font_family)),
             (max_x * 0.95, start_y - 2 * annotation_spacing, text("Iterations: $formatted_iterations", font_color, :right, font_size, font_family)),
             (max_x * 0.95, start_y - 3 * annotation_spacing, text("Error (MSE): $(fvals[itr]) mV", font_color, :right, font_size, font_family))])



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

using Plots

function plotresults2(pr::NamedTuple, res::NamedTuple; savePlot=true, saveLocation="processedData/")
    @unpack converged, statusMessage, fvals, xvals, backtrackVals = res
    n = size(xvals, 1)
    itr = size(xvals, 2)
    x☆ = xvals[:, itr]
    minBacktracks, maxBacktracks = extrema(backtrackVals)

    t_values = pr.df.t
    V_values = pr.df.V
    
    model_vals = [pr.objective(x☆, t, getGradientToo=false) for t in t_values]
    println(length(t_values))
    println(length(model_vals))
    p1 = scatter(t_values, V_values, label="Data", color=:purple, legend=false)
    plot!(p1, t_values, model_vals, label="Model", color=:blue, linewidth=2)

    title_text = "Voltage vs Time\n"
    method_text = "Method: $(pr.alg.method)\n"
    line_search_text = "LineSearch: $(pr.alg.linesearch)\n"
    iterations_text = "Iterations: $(itr)\n"
    error_text = "Error (MSE): $(fvals[itr]) mV"

    p2 = plot(title=title_text, framestyle=:none, axis=false)
    p3 = plot(title=method_text, framestyle=:none, axis=false)
    p4 = plot(title=line_search_text, framestyle=:none, axis=false)
    p5 = plot(title=iterations_text, framestyle=:none, axis=false)
    p6 = plot(title=error_text, framestyle=:none, axis=false)

    p7 = plot(p1, p2, p3, p4, p5, p6, layout=(6,1), size=(600,600), margin=5)

    if savePlot
        filename = saveLocation * string(pr.objective) * "_$(pr.alg.method)_$(@sprintf("%.0e", itr)).png"
        savefig(p7, filename)
    end
end

