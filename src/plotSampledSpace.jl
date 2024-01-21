using Plots

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
