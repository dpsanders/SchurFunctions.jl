using DoubleDouble

convert{T<:BitsFloat}(::Type{Double{T}}, x::Int64) = Double(Float64(x))

convert{T<:BitsFloat}(::Type{Double{T}}, x::Real) = Double(float(x))

Base.round(Int, x::Double) = round(Int, Float64(x))
