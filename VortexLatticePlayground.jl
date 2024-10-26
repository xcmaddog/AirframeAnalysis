using VortexLattice
using Plots
using StaticArrays
include("MyVortexLattice.jl")
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

#=
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

CD, CY, CL = CF # CD = drag, CY = side forces, CL = lift
Cl, Cm, Cn = CM # Cl = , Cm = pitching moment, Cn =

MyVortexLattice.draw_airframe([wing_grid, tail_grid])

MyVortexLattice.print_results(CF, CM, CDiff)

properties = get_surface_properties(system)
write_vtk("playgroundwing", surfaces, properties, symmetric = [true, true])
=#

#=
#define parameters
wing_root_chord = 5
wing_span = 20
wing_angle = 2
spanwise_sections = 10:20

# Set up freestream parameters
alpha = -1.0 * pi / 180 # angle of attack
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
Vinf = 1.0 # reference velocity
fs = Freestream(Vinf, alpha, beta, Omega)

# Initialize reference location and symmetric option
rref = SVector(0.50, 0.0, 0.0) # using StaticArrays for fixed-size vector
symmetric = true

# Initialize vectors to store results for each coefficient
CDs, CYs, CLs = Float64[], Float64[], Float64[]
Cls, Cms, Cns = Float64[], Float64[], Float64[]
CDiffs = Float64[]

# Run analysis for each spanwise section count
for n in spanwise_sections
    # geometry for the right half of the wing
    wing_xyz, wing_reference_area, wing_reference_chord = 
        MyVortexLattice.grid_from_elliptical_edge(wing_root_chord, wing_span, 6, n, wing_angle, 0, 0)

    # reference parameters for each section count
    ref = Reference(wing_reference_area, wing_reference_chord, wing_span, rref, Vinf)

    # construct surface
    wing_grid, wing_surface = grid_to_surface_panels(wing_xyz, mirror = false)
    
    # create vector containing the surface and perform steady state analysis
    surfaces = [wing_surface]
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
    
    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    # retrieve far-field drag
    CDiff = far_field_drag(system)
    
    # Append the results for each coefficient
    push!(CDs, CD)
    push!(CYs, CY)
    push!(CLs, CL)
    push!(Cls, Cl)
    push!(Cms, Cm)
    push!(Cns, Cn)
    push!(CDiffs, CDiff)
end

# At this point, the six vectors contain the results for each spanwise section count
println("CDs: ", CDs)
println("CYs: ", CYs)
println("CLs: ", CLs)
println("Cls: ", Cls)
println("Cms: ", Cms)
println("Cns: ", Cns)
println("CDiffs: ", CDiffs)
=#

#define wing parameters
wing_root_chord = 5
wing_span = 20
wing_angle = 5
spanwise_sections = 60

#define horizontal tail parameters
tail_root_chord = 2
horizontal_tail_span = 5
tail_x_placement = 15
tail_z_placement = 1

#define vertical tail parameters
#initial_vertical_tail_span = 0.01
initial_vertical_tail_span = 5

# Set up freestream parameters
alpha = -1.0 * pi / 180 # angle of attack
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
Vinf = 1.0 # reference velocity
fs = Freestream(Vinf, alpha, beta, Omega)

# Initialize reference location and symmetric option
rref = SVector(0.50, 0.0, 0.0) # using StaticArrays for fixed-size vector
symmetric = false

# geometry for the right half of the wing
wing_xyz, wing_reference_area, wing_reference_chord = 
    MyVortexLattice.grid_from_elliptical_edge(wing_root_chord, wing_span, 6, spanwise_sections, wing_angle, 0, 0, false)

#create the horizontal tail geometry
horizontal_tail_xyz, horizontal_tail_reference_area, horizontal_tail_reference_chord = 
    MyVortexLattice.grid_from_elliptical_edge(tail_root_chord, horizontal_tail_span, 6, spanwise_sections, 0, tail_x_placement, tail_z_placement, false)

# construct wing surface
wing_grid, wing_surface = grid_to_surface_panels(wing_xyz, mirror = true)

#construct the horizontal tail surface
horizontal_tail_grid, horizontal_tail_surface = grid_to_surface_panels(horizontal_tail_xyz, mirror = true)

# reference parameters (determined by the wing)
ref = Reference(wing_reference_area, wing_reference_chord, wing_span, rref, Vinf)

# Initialize vectors to store results for each derivative
dCFs, dCMs = [], []
ratios = []

# Run analysis for each spanwise section count
for n in range(0,2)
    vertical_tail_span = initial_vertical_tail_span + (0.01 * n)

    #geometry for vertical tail
    vertical_tail_xyz, vertical_tail_reference_area, vertical_tail_reference_chord = 
        MyVortexLattice.grid_from_elliptical_edge(tail_root_chord, vertical_tail_span, 6, spanwise_sections, 0, tail_x_placement, tail_z_placement, true)

    #construct vertical tail surface
    vertical_tail_grid, vertical_tail_surface = grid_to_surface_panels(vertical_tail_xyz, mirror = false)

    # create vector containing the surfaces and perform steady state analysis
    surfaces = [wing_surface, horizontal_tail_surface, vertical_tail_surface]
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
    
    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    #retreive stability derivatives
    body_derivatives(system)
    dCF, dCM = stability_derivatives(system)

    
    # Append the results for each coefficient
    push!(dCFs, dCF)
    push!(dCMs, dCM)
    push!(ratios, horizontal_tail_reference_area / vertical_tail_reference_area)

    properties = get_surface_properties(system)
    write_vtk("withVerticalTail" * string(n), surfaces, properties, symmetric = [false, false, false])
end