
using CSV
using DataFrames

func_dir = @__DIR__
root_dir = dirname(dirname(func_dir))
rawData_dir = joinpath(root_dir, "rawData")
# Function to create a spiral matrix
function create_spiral_matrix(size, high_speed, low_speed, thickness, loops)
    # Create a grid of low_speed values
    speed_matrix = fill(low_speed, size, size)
    
    # Calculate the number of points in the spiral
    points = size * size ÷ (10 * loops)
    
    # Create the spiral pattern
    theta = range(0, stop=2pi*loops, length=points)
    r = range(0, stop=size÷2, length=points)
    x_spiral = Int.(round.(size÷2 .+ r .* cos.(theta)))
    y_spiral = Int.(round.(size÷2 .+ r .* sin.(theta)))
    
    # Set high_speed values along the spiral path
    for i in 1:length(x_spiral)
        x, y = x_spiral[i], y_spiral[i]
        # Ensure the indices are within the bounds of the matrix
        if x > 0 && x ≤ size && y > 0 && y ≤ size
            for dx in -thickness÷2:thickness÷2
                for dy in -thickness÷2:thickness÷2
                    xi, yi = x+dx, y+dy
                    if xi > 0 && xi ≤ size && yi > 0 && yi ≤ size
                        speed_matrix[yi, xi] = high_speed
                    end
                end
            end
        end
    end
    
    return speed_matrix
end

# Parameters for the spiral matrix
size = 256
high_speed = 9.0
low_speed = 1.0
thickness = 10 # Thickness of the spiral lines
loops = 5  # Number of loops in the spiral

# Generate the spiral speed matrix
spiral_speed_matrix = create_spiral_matrix(size, high_speed, low_speed, thickness, loops)

filename_csv = joinpath(rawData_dir, "SpeedData_Spiral.csv")
# Save the matrix to a CSV file
# CSV.write("SpeedData_Spiral.csv", DataFrame(spiral_speed_matrix, :auto))
CSV.write(filename_csv, DataFrame(spiral_speed_matrix, :auto), header=false)
