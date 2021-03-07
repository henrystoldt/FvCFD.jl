include("fvCFD.jl")

# ### Delta Wing ###
# # Freestream Conditions
# P = 7669
# T = 180
# U = [ 539, 0, 0 ]
#
# # Choose a mesh
# meshPath = "../test/OFmemesDeltaMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], wallBoundary, [] ]