# Madsen's Wrapper for Vortex Lattice
#=
(a) Write a function that takes in the root chord and span of a wing as well as a number
of sections used to define the wing.
(b) Have the function use the root chord and span numbers to define an ellipse.
(c) Use the number of sections to divide up the wing span.
(d) Have the function return the values of the leading edge locations and chord lengths
at each of the sections so that your wing geometry approximates an ellipse.
(e) Compare the efficiency of wings generated with only a few sections to those generated
with many sections.
=#

using Plots
using VortexLattice

module MyVortexLattice

export ellipicalWingGenerator, grid_from_elliptical_edge

"""
ellipicalWingGenerator(root_chord, span, num_of_sections)
This function calculates the positions of points along the leading edge
    of the right side of an ellipical wing
root_chord defines the maximum length of the chord
span is the span of the entire wing (not just the right half)
num_of_sections is the number of points defined along the leading edge
"""
function ellipicalWingGenerator(root_chord, span, num_of_sections)
    #this function calculates the positions of wing tip locations for
        #an apporoximation of an ellipical wing. Note that the coordinate
        #system used places the orgin in the center of the wing
    a = root_chord / 2
    b = span / 2
    y = range(start = 0, stop = b, length = num_of_sections)
    x = (a ^ 2) .* (1 .- ((y .^ 2) ./ (b ^ 2)))
    x = sqrt.(x)
    chord_lengths = 2 .* x
    return x, y, chord_lengths
end
"""
grid_from_elliptical_edge(root_chord, span, chord_points, span_points)
This function returns a grid of points that represent the right side of a level wing
root_chord defines the maximum length of the chord
span is the span of the entire wing (not just the right half)
chord_points is the number of chordwise points
span_points is the number of spanwise points
This function returns a (3, chord_points, span_points) matrix of points
"""
function grid_from_elliptical_edge(root_chord, span, chord_points, span_points)
    #initialize the grid
    grid = zeros(Float64, 3, chord_points, span_points)
    
    x_leading, y_leading, chord_lengths = ellipicalWingGenerator(root_chord, span, span_points)

    for j in 1:span_points #traverse spanwise
        for i in 1:chord_points #traverse chordwise
            #define x
            grid[1, i, j] = x_leading[j] - (((i-1) / (chord_points-1)) * chord_lengths[j])
            #define y
            grid[2, i, j] = y_leading[j]
        end
    end
    return grid
end

end; #this is the end for the module

#=
import .MyVortexLattice

grid = grid_from_elliptical_edge(5, 20, 5, 15)

scatter3d(grid[1, 1 ,:], grid[2,1,:], grid[3,1,:])

for i in range(2,5)
    scatter3d!(grid[1, i ,:], grid[2, i,:], grid[3, i,:])
end

#plot!(display)
#savefig("myPlot.png")
=#