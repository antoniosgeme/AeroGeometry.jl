using AeroGeometry
using TestItemRunner
using Test

@run_package_tests verbose=true

include("test_airfoils.jl")

include("test_wings.jl")

include("test_readme.jl")