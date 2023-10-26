
using DataFrames
using LaTeXStrings
using Parameters
using Plots
using Printf

# include("utilities.jl")
include("objfuns/objective.jl")

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

function showresults(res::NamedTuple;
    log::Bool=true,
    log_path::String="./logging/")

    @unpack converged, statusMessage, fvals, αvals, backtrackVals, xvals, M, fevals, gevals, cause, pr = res

    dataFitting = isempty(pr.p.data) ? false : true

    result_txt = log_path*"results_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    if isfile(result_txt)
        rm(result_txt) # remove results log file if already present
    end
    log, log_path = log, result_txt  # or whatever default values you have

    function myprintln1(v, message)
        myprintln(v, message, log=log, log_path=result_txt)
    end
    # println(res.pr)
    
    nnztol = 1e-8 # variable having a value below this will NOT be counted in the list of non-zero variables. For printing purposes only. 

    v = true
    myprintln1(v, "****************************")
    myprintln1(true, "Solver run has concluded.")
    myprintln1(true, "Function solved for: $(pr.objective)()")
    if dataFitting
        myprintln1(true, "This was a data fitting problem.")
    else
        myprintln1(true, "This was a function minimization problem.")
    end
    myprintln1(true, "Cause(s) for stopping:")
    causeForStopping = res.cause
    for msg ∈ causeForStopping
        myprintln1(true, msg)
    end
    if converged == true
        myprintln1(true, statusMessage)
    elseif converged == false
        myprintln1(true, statusMessage)
    else
        @error "bad condition"
    end

    methodMsg = "Method Used: " * pr.alg.method
    lsMsg = "Linesearch method used: " * pr.alg.linesearch
    myprintln1(v, "****************************")
    myprintln1(v, methodMsg)
    myprintln1(v, lsMsg)

    
    n = size(xvals, 1)
    itr = size(xvals, 2)
    x☆ = xvals[:, itr]
    minBacktracks, maxBacktracks = extrema(backtrackVals)

    myprintln1(v, "****************************")
    
    if converged == true
        fvalPrefixMSE = "Optimal MSE value = "
        fvalPrefixSSE = "Optimal SSE value = "
        fvalPrefixMinimum = "Optimal function value = "
        xvalMessage = "Optimal Nonzero Variables:"
    elseif converged == false
        fvalPrefixMSE = "Best MSE value = "
        fvalPrefixSSE = "Best SSE value = "
        fvalPrefixMinimum = "Best function value = "
        xvalMessage = "Best Known Nonzero Variables:"
    else
        @error "bad condition"
    end
    fvalMessageMSE = fvalPrefixMSE*"$(fvals[itr])"
    fvalMessageSSE = fvalPrefixSSE*"$(fvals[itr]*M)"
    fvalMessageMinimum = fvalPrefixMinimum*"$(fvals[itr])"

    if dataFitting
        myprintln1(v, fvalMessageMSE)
        myprintln1(v, fvalMessageSSE)
    else
        myprintln1(v, fvalMessageMinimum)
    end

    myprintln1(v, "***************************")
    myprintln1(v, xvalMessage)
    for k in 1:n
        if abs(x☆[k]) > nnztol     
            myprintln1(v, "x☆[$(k)] = $(x☆[k])")
        end
    end
    myprintln1(v, "***************************")
    myprintln1(v, "Number of fevals = $(fevals)")
    myprintln1(v, "Number of gevals = $(gevals)")
    myprintln1(v, "Number of evals = $(fevals+gevals*n)")

    myprintln1(v, "*"^50)
    myprintln1(v, "Tolerance critera used:")
    myprintln1(v, "dftol = $(res.pr.alg.dftol)")

    # only if xopt given, below xopt is for rosenbrock function: a vector of ones
    # diffx = sum( abs.( res.xvals[:, end] - ones(Float64, length(res.xvals[:, end])) ) )
end


using LaTeXStrings
using Plots

function plotresults(pr::NamedTuple, res::NamedTuple; savePlot::Bool=true, saveLocation::String="processedData/")
    @unpack converged, statusMessage, fvals, xvals, backtrackVals = res
    n = size(xvals, 1)
    itr = size(xvals, 2)
    x☆ = xvals[:, itr]

    # Convert time to milliseconds
    time_ms = pr.df.x .* 1000
    
    # Calculate dampedSHM values using optimal parameters
    predicted_values = [pr.objective(x☆, t, getGradientToo=false) for t in pr.df.x]

    # Plot
    gr()  # Ensure you're using the GR backend for Plots.jl which supports LaTeX rendering
    theme(:dark)  # Setting a dark theme

    formatted_iterations = @sprintf("%.0e", itr)

    title_content = "Voltage vs Time\nMethod: $(pr.alg.method)\nLineSearch: $(pr.alg.linesearch)\nIterations: $formatted_iterations\nError (MSE): $(fvals[itr]) mV"

    p = scatter(time_ms, pr.df.y, 
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
max_x = maximum(pr.df.x)
max_y = minimum(pr.df.y)
annotation_spacing = (maximum(pr.df.y) - max_y) * 0.05

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

