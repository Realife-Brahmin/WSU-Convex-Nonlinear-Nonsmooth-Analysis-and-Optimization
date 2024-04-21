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