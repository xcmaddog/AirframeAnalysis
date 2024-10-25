using VortexLattice
using Plots
include("MyVortexLattice.jl")
#using .MyVortexLattice
import .MyVortexLattice

"""
testPoints()
This function runs a simple test on MyVortexLattice.grid_from_elliptical_edge
    and outputs a graph representing the grid created
"""
function testPoints()
    grid, -, - = MyVortexLattice.grid_from_elliptical_edge(5, 20, 6, 15, 0)

    scatter(grid[1, 1 ,:], grid[2,1,:])

    for i in range(2,6)
        scatter!(grid[1, i ,:], grid[2, i,:])
    end
    scatter!(show = true)
end

#testPoints()

#define parameters
wing_root_chord = 5
wing_span = 20
wing_angle = 2

tail_root_chord = 2
tail_span = 5
tail_x_placement = 15
tail_z_placement = 3

#root_chord, span, chord_points, span_points, rotation, chordwise_translation, vertical_translation
# geometry (right half of the wing)
wing_xyz, wing_refrence_area, wing_refrence_chord = 
    MyVortexLattice.grid_from_elliptical_edge(wing_root_chord, wing_span, 6, 15, wing_angle, 0, 0)
 
# geometry (right half of tail)
tail_xyz, tail_reference_area, tail_reference_chord = 
    MyVortexLattice.grid_from_elliptical_edge(tail_root_chord, tail_span, 6, 10, 0, tail_x_placement, tail_z_placement)

# reference parameters (These are based on the wing only, excluding the tail)
#Sref = 30.0 # reference area
#cref = 2.0  # reference chord
#bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(wing_refrence_area, wing_refrence_chord, wing_span, rref, Vinf)

# freestream parameters
alpha = -1.0*pi/180 # angle of attack (is negative because here we are defining the angle of the freestream, not the wing)
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
#grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
#    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)
wing_grid, wing_surface = grid_to_surface_panels(wing_xyz, mirror = false)
tail_grid, tail_surface = grid_to_surface_panels(tail_xyz, mirror = false)

# create vector containing all surfaces
surfaces = [wing_surface, tail_surface]
#surfaces = [wing_surface]
#surfaces = [tail_surface]

# we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
symmetric = true

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system)

#CD, CY, CL = CF # CD = drag, CY = side forces, CL = lift
#Cl, Cm, Cn = CM # Cl = , Cm = pitching moment, Cn =


function printResults(CF, CM, CDiff)
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

function draw_airframe(grids)
    scatter3d(aspect_ratio=:equal)
    for grid in grids
        for i in range(1, size(grid)[2])
            for j in range(1, size(grid)[3])
                scatter3d!([grid[1, i, j]], [grid[2, i, j]], [grid[3, i, j]])
            end
        end
    end
    scatter3d!(show = true)
end

#draw_airframe([wing_grid, tail_grid])

printResults(CF, CM, CDiff)

properties = get_surface_properties(system)
write_vtk("playgroundwing", surfaces, properties, symmetric = [true, true])