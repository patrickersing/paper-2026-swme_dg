using Trixi
using TrixiShallowWater
using Static: True, False
using StaticArrays: SVector, MVector, @SMatrix, SMatrix, MArray, SArray
using MuladdMacro: @muladd
using LinearAlgebra: diagm

include("equations/equations.jl")
include("equations/moment_matrices.jl")
