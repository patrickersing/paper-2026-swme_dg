using Trixi
using TrixiShallowWater
using Static: True, False
using StaticArrays: SVector, MVector, @SMatrix, MArray, SArray
using MuladdMacro: @muladd

include("equations/equations.jl")
include("equations/moment_matrices.jl")
