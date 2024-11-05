#=
Author: Madsen Evans
FLOW Lab BYU
Date: 11/5/2024

This is a module that has some functions useful for creating wings to use
    with the VortexLattice package
=#

module MyVortexLattice

using Rotations
using Plots

export ellipical_wing_generator, grid_from_elliptical_edge, print_results

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
function ellipical_wing_generator(root_chord, span, num_of_sections, chordwise_translation)
    # a and b define the two chords of the ellipse
    a = root_chord / 2
    b = span / 2
    #get the points for the leading edge of the elliptical wing
    y = range(start = 0, stop = b, length = num_of_sections + 1)
    x = (a ^ 2) .* (1 .- ((y .^ 2) ./ (b ^ 2)))
    x = sqrt.(x)
    #store the lengths of the chords for each y point
    chord_lengths = 2 .* x 
    #put the center of the leading edge of the wing at the origin, and then shift it by the chordwise_translation
    x = x .+ (chordwise_translation + a)
    #approximate the area of the wing for reference
    area = 0
    for i in 1:(length(chord_lengths) - 1)
        area = area + (chord_lengths[i] * (y[i+1] - y[i]))
    end
    #calculate a rough average chord for reference
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

Note that the spanwise section representing the tip of the wing is omitted because the chord length is zero
"""
function grid_from_elliptical_edge(root_chord, span, chord_points, span_points, rotation, chordwise_translation, vertical_translation, vertical)
    span_points = span_points + 1
    #initialize the grid
    grid = zeros(Float64, 3, chord_points, span_points)
    
    x_leading, y_leading, chord_lengths, area, average_chord = ellipical_wing_generator(root_chord, span, span_points, chordwise_translation)

    for j in 1:span_points #traverse spanwise
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

    if vertical 
        rot = Rotations.RotX(-pi/2)
        for j in 1:size(grid)[3] #traverse spanwise
            for i in 1:size(grid)[2] #traverse chordwise
                point = [grid[1,i,j] grid[2,i,j] grid[3,i,j]] #pull the point into a single vector
                point = point * rot #rotate the point
                grid[1,i,j] = point[1] #put the point back ino the grid
                grid[2,i,j] = point[2]
                grid[3,i,j] = point[3]
            end
        end
        #move it to where it should be
        for j in 1:size(grid)[3] #traverse spanwise
            for i in 1:size(grid)[2] #traverse chordwise
                grid[3,i,j] = grid[3,i,j] + vertical_translation
                grid[2,i,j] = 0
            end
        end

    else
        area = area * 2 #we are assuming that vertical wings aren't mirrored
    end

    return grid, area, average_chord
end

"""
printResults(CF, CM, CDiff)

This is a simple function that takes the resulting flight coefficients and prints them out
"""
function print_results(CF, CM, CDiff)
    println("Coefficient of drag: " * string(round(CF[1], digits = 6)))
    println("Coefficient of spanwise forces: " * string(round(CF[2], digits = 6)))
    println("Coefficient of lift: " * string(round(CF[3], digits = 6)))
    println()
    println("Coefficient of roll or yaw: " * string(round(CM[1], digits = 6)))
    println("Coefficient of pitching moment: " * string(round(CM[2], digits = 6)))
    println("Coefficient of roll or yaw: " * string(round(CM[3], digits = 6)))
    println()
    println("Coefficient of drag from far field analysis: " * string(round(CDiff, digits = 6)))
end

"""
draw_airframe(grids)

This function plots the wings as a 3d scatter plot

grids is a vector of (3,i,j) matricies where i is chordwise and j is spanwise
"""
function draw_airframe(grids)
    scatter3d(aspect_ratio=:equal) #prep the plot
    for grid in grids #for each wing
        for i in range(1, size(grid)[2]) #traverse chordwise
            for j in range(1, size(grid)[3]) #traverse spanwise
                scatter3d!([grid[1, i, j]], [grid[2, i, j]], [grid[3, i, j]]) #plot the point
            end
        end
    end
    scatter3d!(show = true)#display the scatter plot
end

end; #this is the end for the module
