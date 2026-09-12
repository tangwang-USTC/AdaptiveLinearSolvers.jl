abstract type AbstractLinearOperator end

"""Linear operator wrapper that counts applications and enforces an optional hard limit."""
mutable struct CountingOperator{T, A} <: AbstractLinearOperator
    source::A
    applications::Int
    limit::Int
end

function CountingOperator(source; limit::Integer=0)
    limit >= 0 || throw(ArgumentError("operator application limit must be nonnegative"))
    return CountingOperator{eltype(source), typeof(source)}(source, 0, Int(limit))
end

Base.size(operator::CountingOperator) = size(operator.source)
Base.size(operator::CountingOperator, dimension::Integer) = size(operator.source, dimension)
Base.eltype(::Type{CountingOperator{T, A}}) where {T, A} = T
Base.eltype(operator::CountingOperator{T}) where {T} = T

function LinearAlgebra.mul!(y::AbstractVector, operator::CountingOperator, x::AbstractVector)
    operator.limit > 0 && operator.applications >= operator.limit &&
        throw(OperatorApplicationBudgetExceeded(operator.limit))
    operator.applications += 1
    return LinearAlgebra.mul!(y, operator.source, x)
end

function Base.:*(operator::CountingOperator{T}, x::AbstractVector) where {T}
    y = similar(x, promote_type(T, eltype(x)), size(operator, 1))
    return LinearAlgebra.mul!(y, operator, x)
end

"""Matrix-free linear operator defined by an in-place action `apply!(y, x)` for `y = A * x`."""
struct MatrixFreeOperator{T, F} <: AbstractLinearOperator
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
Base.size(operator::MatrixFreeOperator, dimension::Integer) =
    dimension == 1 ? operator.rows : dimension == 2 ? operator.cols : 1
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

"""Contiguous block partition of a flat vector used by `BlockOperator`."""
struct BlockLayout
    sizes::Vector{Int}
    offsets::Vector{Int}
end

function BlockLayout(sizes::AbstractVector{<:Integer})
    block_sizes = Int[sizes...]
    isempty(block_sizes) && throw(ArgumentError("a block layout requires at least one block"))
    all(size -> size > 0, block_sizes) || throw(ArgumentError("block sizes must be positive"))
    offsets = cumsum(vcat(0, block_sizes))
    return BlockLayout(block_sizes, offsets)
end

blockrange(layout::BlockLayout, block::Integer) =
    (layout.offsets[block] + 1):layout.offsets[block + 1]

"""Block-composed matrix-free operator acting on flat vectors with a declared block layout."""
struct BlockOperator{T, B} <: AbstractLinearOperator
    layout::BlockLayout
    blocks::B
end

function BlockOperator(layout::BlockLayout, blocks::AbstractMatrix; T::Type=Float64)
    number_of_blocks = length(layout.sizes)
    size(blocks) == (number_of_blocks, number_of_blocks) ||
        throw(DimensionMismatch("block array shape must match the block layout"))
    for row in 1:number_of_blocks, column in 1:number_of_blocks
        block = blocks[row, column]
        block === nothing && continue
        expected_size = (layout.sizes[row], layout.sizes[column])
        size(block) == expected_size ||
            throw(DimensionMismatch("block ($row, $column) has size $(size(block)); expected $expected_size"))
    end
    return BlockOperator{T, typeof(blocks)}(layout, blocks)
end

Base.size(operator::BlockOperator) = (operator.layout.offsets[end], operator.layout.offsets[end])
Base.size(operator::BlockOperator, dimension::Integer) =
    dimension == 1 || dimension == 2 ? operator.layout.offsets[end] : 1
Base.eltype(::Type{BlockOperator{T, B}}) where {T, B} = T
Base.eltype(operator::BlockOperator{T}) where {T} = T

function LinearAlgebra.mul!(y::AbstractVector, operator::BlockOperator, x::AbstractVector)
    total_size = operator.layout.offsets[end]
    length(y) == total_size || throw(DimensionMismatch("output length does not match block layout"))
    length(x) == total_size || throw(DimensionMismatch("input length does not match block layout"))
    result_type = promote_type(eltype(y), eltype(x))
    number_of_blocks = length(operator.layout.sizes)
    for row in 1:number_of_blocks
        y_block = view(y, blockrange(operator.layout, row))
        fill!(y_block, zero(eltype(y_block)))
        for column in 1:number_of_blocks
            block = operator.blocks[row, column]
            block === nothing && continue
            x_block = view(x, blockrange(operator.layout, column))
            temporary = similar(y_block, result_type, length(y_block))
            LinearAlgebra.mul!(temporary, block, x_block)
            y_block .+= temporary
        end
    end
    return y
end

function Base.:*(operator::BlockOperator{T}, x::AbstractVector) where {T}
    result_type = promote_type(T, eltype(x))
    y = similar(x, result_type, size(operator, 1))
    return LinearAlgebra.mul!(y, operator, x)
end
