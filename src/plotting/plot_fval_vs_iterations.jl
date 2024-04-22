"""
    plot_fval_vs_iterations(res; savePlot=true)

Plot the progression of function values against iterations for a given optimization result or multiple runs.

# Arguments
- `res`: The result or a vector of results from optimization processes. Each result must contain `pr`, `fvals`, and `xvals` as part of its structure.

# Keyword Arguments
- `savePlot::Bool=true`: Flag to determine whether to save the plot to a file. If true, the plot is saved under the processed data directory.

# Returns
- `p1`: A plot object representing the function value progression over iterations.

# Description
This function visualizes the optimization performance by plotting function values against iteration numbers. It supports handling either a single result or a collection of results from multiple optimization runs. When given multiple results, it recursively calls itself to plot each individual result. 

The function extracts necessary data such as function values and iteration counts from the optimization result and uses this data to create a plot. For single runs, it plots the function value over each iteration, highlighting the minimum function value achieved and the corresponding iteration.

# Notes
- The plot includes a title that provides information about the minimum function value reached, the optimization method used, and other relevant details.
- The plots are saved in a structured directory format using the base directory `dirname(dirname(@__DIR__))` appended with `/processedData/`.
- It is important to ensure that each result in `res` contains all necessary data (`pr`, `fvals`, `xvals`) for the function to process and plot the information correctly.
- The function employs the `Plots` library with the `GR` backend specified for plotting.
- Plot aesthetics such as fonts, labels, and line properties are customized to enhance readability and presentation quality.

# Example
```julia
# Assuming `results` is a single optimization result or a vector of results
plot_fval_vs_iterations(results, savePlot=true)

# This will generate a plot for each optimization run in `results`, displaying function values against iterations and saving the plots if `savePlot` is true.
```
"""
function plot_fval_vs_iterations(res;
    savePlot=true)

        if typeof(res) <: Vector
            numRuns = length(res)
            for runNum ∈ 1:numRuns
                res1 = res[runNum]
                plot_fval_vs_iterations(res1)
            end

            return

        else

            res1 = res

        end

        pr = res1.pr
        functionName = pr.objectiveString
        theme(:dao)
        fvals = res1.fvals
        itr = length(fvals)
        if itr == 0 # x0 is x☆
            x0 = pr.x0
            pDict = pr.p
            f0 = pr.objective(x0, pDict, getGradientToo=false)
            x☆ = x0
            f☆ = f0
            x = x☆
            f = f☆
        else
            x☆ = res1.xvals[:, itr]
            f☆ = res1.fvals[itr]
            f = f☆
            x = x☆
        end

        n = length(x)

    
        gr()
        p1 = Plots.plot(1:itr, fvals,
            title = "Function Value vs Iteration: $(L"f_{min} =") $(round(f, digits=5)) \n using $(pr.alg[:method])",
            titlefont = font(12,"Computer Modern"),
            guidefont = font(15,"Computer Modern"),
            label = L"f",
            ylabel = L"f",
            xlabel = L"Iteration",
            linewidth = 3,
            markersize = 2,
            alpha = 1.0,
        )
    
        folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
        filename = folderName*functionName*"_fval_vs_itr_"*pr.alg[:method]*"_"*string(n)*".png"
    
    
        if savePlot
            
            if isfile(filename)
                rm(filename)
            end
        
            Plots.savefig(p1, filename)
        end
    
        return p1
end