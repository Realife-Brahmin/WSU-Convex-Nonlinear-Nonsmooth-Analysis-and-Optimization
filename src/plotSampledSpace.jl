using Plots

"""
    plotSampledSpace(ss, method; view::Bool=true, savePlot::Bool=true)

Plot a sampled space of data points in 2D.

This function generates a scatter plot to visualize a sampled space of data points in two dimensions (R^2). It provides options to control whether the plot is displayed and whether it is saved as an image.

# Arguments:
- `ss::Matrix{Float64}`: A 2D matrix containing the sampled data points. Each column of `ss` represents a data point with two dimensions (x and y coordinates).
- `method::String`: A string indicating the method used for sampling the data.
- `view::Bool=true`: (Optional) If `true`, the plot is displayed. If `false`, the plot is not displayed. Default is `true`.
- `savePlot::Bool=true`: (Optional) If `true`, the plot is saved as an image. If `false`, the plot is not saved. Default is `true`.

# Example:
```julia
# Create a sampled space matrix
ss = rand(2, 100)

# Plot the sampled space using the "Halton" method, display it, and save it as an image
plotSampledSpace(ss, "Halton", view=true, savePlot=true)
"""
function plotSampledSpace(ss, method;
        view::Bool=true,
        savePlot::Bool=true)

        n, p = size(ss)
        if n != 2
            @error "Can only plot planar plot for R^2!"
        end
        Plots.theme(:dao)
        p1 = scatter(ss[1, :], ss[2, :],
                aspect_ratio=:equal,
                title="Sampled Space using $(method) method",
                xlims=[0.0, 1.0],
                ylims=[0.0, 1.0],
                label=:none);
        
        if view        
            display(p1)
        end

        if savePlot
            ext = ".png"
            filename = method*"_n_"*string(n)*"_p_"*string(p)*ext
            fullAdd = joinpath("processedData", filename)
            savefig(p1, fullAdd)
        end
end;
