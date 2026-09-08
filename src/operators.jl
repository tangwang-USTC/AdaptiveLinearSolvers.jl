"""Matrix-free linear operator defined by an in-place action `apply!(y, x)` for `y = A * x`."""
struct MatrixFreeOperator{T, F}
    rows::Int
    cols::Int
    apply!::F
end

function MatrixFreeOperator(rows::Integer, cols::Integer, apply!; T::Type=Float64)
    rows >= 0 || throw(ArgumentError("rows must be nonnegative"))
    cols >= 0 || throw(ArgumentError("cols must be nonnegative"))
    return MatrixFreeOperator{T, typeof(apply!)}(Int(rows), Int(cols), apply!)
end

Base.size(operator::MatrixFreeOperator) = (operator.rows, operator.cols)
Base.eltype(::Type{MatrixFreeOperator{T, F}}) where {T, F} = T
Base.eltype(operator::MatrixFreeOperator{T}) where {T} = T

function LinearAlgebra.mul!(y::AbstractVector, operator::MatrixFreeOperator, x::AbstractVector)
    length(y) == operator.rows || throw(DimensionMismatch("output length does not match operator rows"))
    length(x) == operator.cols || throw(DimensionMismatch("input length does not match operator columns"))
    operator.apply!(y, x)
    return y
end

function Base.:*(operator::MatrixFreeOperator{T}, x::AbstractVector) where {T}
    result_type = promote_type(T, eltype(x))
    y = similar(x, result_type, operator.rows)
    return LinearAlgebra.mul!(y, operator, x)
end
