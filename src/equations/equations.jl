
abstract type AbstractMomentEquations{NDIMS, NVARS, NMOMENTS} <:
              Trixi.AbstractEquations{NDIMS, NVARS} end

"""
    eachmoment(equations::AbstractMomentEquations)
    
Return an iterator over the indices that specify the location in relevant data structures
for the moments in `AbstractMomentEquations`. 
"""
@inline function eachmoment(equations::AbstractMomentEquations)
    Base.OneTo(nmoments(equations))
end

"""
    nmoments(equations::AbstractMomentEquations)

Retrieve the number of moments from an equation instance of the `AbstractMomentEquations`.
"""
@inline function nmoments(::AbstractMomentEquations{NDIMS, NVARS, NMOMENTS}) where {NDIMS,
                                                                                    NVARS,
                                                                                    NMOMENTS
                                                                                    }
    NMOMENTS
end

include("shallow_water_linearized_moments_1d.jl")
include("shallow_water_moments_1d.jl")
