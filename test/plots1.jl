using Plots

include("rawDataTest/medianDevSalaries.jl"); # contains several data vectors, which can direcly be used for plotting

theme(:dracula)
p1 = plot(dev_x1, dev_y1, dpi = 7000, 
    xlabel = "Developer Age [years]",
    ylabel = "Median Developer Salaries [\$]",
    title = "Developer Salaries",
    label = "Gen Developer Salaries",
    color = :chocolate1,
    # color = :slateblue,
    # color = :gray2,
    # color = :maroon3,
    alpha = 1.0,
    linestyle = :solid,
    linewidth = 3.0)

p2 = plot(p1, dev_x1, py_dev_y1, dpi = 7000, 
label = "Python Developer Salaries",
# color = :orchid,
color = :teal,
# color = :gray2,
# color = :maroon3,
alpha = 0.9,
linestyle = :solid,
linewidth = 4.5) 

p3 = plot(p2, dev_x1, js_dev_y1, dpi = 7000, 
label = "Javascript Developer Salaries",
# color = :orchid,
# color = :teal,
color = :antiquewhite1,
# color = :gray2,
# color = :maroon3,
linestyle = :dash,
alpha = 0.8,
linewidth = 2.5)