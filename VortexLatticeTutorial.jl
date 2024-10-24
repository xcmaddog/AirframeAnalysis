using VortexLattice
include("MyVortexLattice.jl")
using .MyVortexLattice
import .MyVortexLattice

# geometry (right half of the wing)
xle = [0.0, 0.4] # leading edge x-position --backwards
yle = [0.0, 7.5] # leading edge y-position --to the right 
zle = [0.0, 0.0] # leading edge z-position -- up (is positive for right hand coordinate system)
chord = [2.2, 1.8] # chord length
theta = [2.0*pi/180, 2.0*pi/180] # twist (in radians)
phi = [0.0, 0.0] # section rotation about the x-axis
fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c))

#xle, yle, chord = ellipicalWingGenerator(4, 15, 12)

# discretization parameters
ns = 12 # number of spanwise panels
nc = 6  # number of chordwise panels
spacing_s = Sine() # spanwise discretization scheme
spacing_c = Uniform() # chordwise discretization scheme

# reference parameters
Sref = 30.0 # reference area
cref = 2.0  # reference chord
bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 1.0*pi/180 # angle of attack
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces = [surface]

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