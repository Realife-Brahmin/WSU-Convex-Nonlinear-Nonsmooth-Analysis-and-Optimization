function scatter_voltage_vs_time(df)
    gr()  # Ensure you're using the GR backend for Plots.jl which supports LaTeX rendering

    theme(:dark)  # Setting a dark theme
    x = df[:, 1]
    y = df[:, 2]
    # Convert time to milliseconds
    time_ms = x .* 1000

    # Create the scatter plot
    scatter(time_ms, y, 
        label=L"Voltage $(V)$",  # Using LaTeXStrings with \text for regular text
        xlabel=L"Time $(ms)$",  # Using LaTeXStrings with \text for regular text
        ylabel=L"Voltage $(mV)$",  # Using LaTeXStrings with \text for regular text
        title="Voltage vs Time",
        markersize=4,
        color=:purple,
        alpha=0.8
    )
end