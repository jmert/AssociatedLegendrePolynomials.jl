import Base.Broadcast: Broadcasted, AbstractArrayStyle, dotview, materialize

"""
    mutable struct Scalar{T} <: AbstractArray{T,0}

A 0-dimensional (scalar) object, much like a `RefValue{T}`, but one which is an
abstract array and assignable under broadcasting operations.
"""
mutable struct Scalar{T} <: AbstractArray{T,0}
    x::T
    Scalar{T}()  where {T} = new()
    Scalar{T}(x) where {T} = new(x)
end
Scalar(s::T) where {T} = Scalar{T}(s)
Base.isassigned(s::Scalar) = isdefined(s, :x)

# AbstractArray interfaces
Base.size(s::Scalar) = ()
@propagate_inbounds Base.getindex(s::Scalar) = s.x
@propagate_inbounds Base.getindex(s::Scalar, ::CartesianIndex{0}) = s.x
@propagate_inbounds Base.setindex!(s::Scalar, x) = (s.x = x; s)
@propagate_inbounds Base.setindex!(s::Scalar, x, ::CartesianIndex{0}) = (s.x = x; s)
Base.similar(s::Scalar{T}, ::Tuple{}) where {T} = Scalar{T}()

# Iteration on a scalar
Base.iterate(s::Scalar) = (s.x, nothing)
Base.iterate(s::Scalar, ::Nothing) = nothing

# Broadcasting extensions
@propagate_inbounds dotview(A::Scalar, ::CartesianIndex{0}) = A
@propagate_inbounds dotview(A::Scalar, ::CartesianIndices{0,Tuple{}}) = A
@propagate_inbounds Base.copyto!(dest::Scalar, bc::Broadcasted{<:AbstractArrayStyle{0}}) =
    dest[] = materialize(bc)
