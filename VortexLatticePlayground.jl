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
    grid = MyVortexLattice.grid_from_elliptical_edge(5, 20, 5, 15)

    scatter(grid[1, 1 ,:], grid[2,1,:])

    for i in range(2,5)
        scatter!(grid[1, i ,:], grid[2, i,:])
    end
    scatter!(show = true)
end

#testPoints()


# geometry (right half of the wing)
xyzw = MyVortexLattice.grid_from_elliptical_edge(5, 20, 6, 15)

# reference parameters
Sref = 30.0 # reference area
cref = 2.0  # reference chord
bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = -4.0*pi/180 # angle of attack (is negative because here we are defining the angle of the freestream, not the wing)
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
#grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
#    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)
grid, wing_surface = grid_to_surface_panels(xyz, mirror = false)

# create vector containing all surfaces
surfaces = [wing_surface]

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

println(CF)
println()
println(CM)
println()
println(CDiff)
