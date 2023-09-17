using Plots
using DataFrames
using LaTeXStrings

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