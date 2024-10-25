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

#using Plots
#using VortexLattice

module MyVortexLattice

using Rotations

export ellipicalWingGenerator, grid_from_elliptical_edge

"""
ellipicalWingGenerator(root_chord, span, num_of_sections)

This function calculates the positions of points along the leading edge
    of the right side of an ellipical wing

root_chord defines the maximum length of the chord
span is the span of the entire wing (not just the right half)
num_of_sections is the number of points defined along the leading edge
chordwise_translation is how far back to put the wing (positive is away from leading edge)

Returns the x and y positions of points along the leading edge, the chord
    length for each section, the area of the entire wing, and the average
    chord length.

Note that the orgin is at the center of the leading edge of the wing before applying the translation.
"""
function ellipicalWingGenerator(root_chord, span, num_of_sections, chordwise_translation)
    a = root_chord / 2
    b = span / 2
    y = range(start = 0, stop = b, length = num_of_sections)
    x = (a ^ 2) .* (1 .- ((y .^ 2) ./ (b ^ 2)))
    x = sqrt.(x)
    chord_lengths = 2 .* x 
    x = x .+ (chordwise_translation + a)
    area = pi * a * b
    average_chord = sum(chord_lengths) / length(chord_lengths)
    return x, y, chord_lengths, area, average_chord
end
"""
grid_from_elliptical_edge(root_chord, span, chord_points, span_points, rotation, chordwise_translation, vertical_translation)

This function returns a grid of points that represent the right side of a level wing

root_chord defines the maximum length of the chord
span is the span of the entire wing (not just the right half)
chord_points is the number of chordwise points
span_points is the number of spanwise points
rotation is an angle in degrees and defines at what angle the wing will be (with reference to the x(chordwise) axis)
chordwise_translation is how far back the wing is moved back
vertical_translation is how far up the wing is moved

This function returns a (3, chord_points, span_points) matrix of points, the
    area of the wing before it was discretized, and an approximation of the 
    average chord length
"""
function grid_from_elliptical_edge(root_chord, span, chord_points, span_points, rotation, chordwise_translation, vertical_translation)
    #initialize the grid
    grid = zeros(Float64, 3, chord_points, span_points-1)
    
    x_leading, y_leading, chord_lengths, area, average_chord = ellipicalWingGenerator(root_chord, span, span_points, chordwise_translation)

    for j in 1:span_points-1 #traverse spanwise
        for i in 1:chord_points #traverse chordwise
            #define x
            grid[1, i, j] = x_leading[j] - (((i-1) / (chord_points-1)) * chord_lengths[j])
            #define y
            grid[2, i, j] = y_leading[j]
        end
    end

    #rotate the wing (about the y/spanwise axis)
    #define the rotation
    rot = Rotations.RotY(-rotation * (pi/180))
    for j in 1:size(grid)[3] #traverse spanwise
        for i in 1:size(grid)[2] #traverse chordwise
            point = [grid[1,i,j] grid[2,i,j] grid[3,i,j]] #pull the point into a single vector
            point = point * rot #rotate the point
            grid[1,i,j] = point[1] #put the point back ino the grid
            grid[2,i,j] = point[2]
            grid[3,i,j] = point[3]
        end
    end

    #translate it back up so that the front of the leading edge is at z = 0 again and translate it any further that it might need
    #get the distance the wing "fell" when rotated
    fall = -grid[3,1,1]

    for j in 1:size(grid)[3] #traverse spanwise
        for i in 1:size(grid)[2] #traverse chordwise
            grid[3,i,j] = grid[3,i,j] + vertical_translation + fall
        end
    end

    return grid, area, average_chord
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