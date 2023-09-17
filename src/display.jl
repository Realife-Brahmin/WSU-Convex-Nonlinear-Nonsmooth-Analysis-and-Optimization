using Plots
using DataFrames
using LaTeXStrings

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

    nnztol = 1e-8
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