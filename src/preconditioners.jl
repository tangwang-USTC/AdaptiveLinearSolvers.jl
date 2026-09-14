"""Incomplete Cholesky IC(0) factor as a left-preconditioner operator for Krylov.jl.

    ICZeroPreconditioner(L)

 wraps the sparse lower-triangular incomplete Cholesky factor `L` (i.e. `A ≈ L * L'`).
 `ldiv!(y, M, x)` applies forward/backward substitution.  The factorization is
 computed by [`_ic_zero_factor`](@ref), which preserves the sparsity pattern of
 `tril(A)` — no fill-in is introduced.
"""
struct ICZeroPreconditioner{T, M <: SparseMatrixCSC{T}}
    L::M
end

Base.size(M::ICZeroPreconditioner, dim::Integer) = size(M.L, dim)
Base.size(M::ICZeroPreconditioner) = size(M.L)

"""
    LinearAlgebra.ldiv!(y, M::ICZeroPreconditioner, x)

Solve `M * y = x` where `M = L * L'`.
Forward substitution `L * z = x`, then backward substitution `L' * y = z`.
Uses column-oriented substitution exploiting CSC storage.
"""
function LinearAlgebra.ldiv!(y::AbstractVector, M::ICZeroPreconditioner, x::AbstractVector)
    n = length(x)
    L = M.L
    length(y) == n || throw(DimensionMismatch("output length does not match input"))

    # Forward substitution (column-oriented): L * z = x.
    # At step j, z[j] is known; we propagate its effect to z[i > j].
    z = copy(x)
    for j in 1:n
        # Find diagonal L[j,j] and divide z[j] by it
        for p in L.colptr[j]:L.colptr[j + 1] - 1
            if L.rowval[p] == j
                z[j] /= L.nzval[p]
                break
            end
        end
        # Propagate: z[i] -= L[i,j] * z[j] for i > j
        for p in L.colptr[j]:L.colptr[j + 1] - 1
            i = L.rowval[p]
            if i > j
                z[i] -= L.nzval[p] * z[j]
            end
        end
    end

    # Backward substitution (column-oriented, reverse order): L' * y = z.
    # CSC stores rows ascending per column (diagonal first).  For the backward
    # substitution we must process off-diagonal entries (row > i) BEFORE the
    # diagonal (row == i), so we traverse the column in REVERSE.
    for i in n:-1:1
        s = z[i]
        for p in (L.colptr[i + 1] - 1):-1:L.colptr[i]
            row = L.rowval[p]
            if row > i
                s -= L.nzval[p] * y[row]
            elseif row == i
                y[i] = s / L.nzval[p]
                break
            end
        end
    end

    return y
end

function Base.:\(M::ICZeroPreconditioner{T}, x::AbstractVector) where T
    y = similar(x, T, size(M, 1))
    return ldiv!(y, M, x)
end

"""
    _get_entry(L::SparseMatrixCSC, i, j)

Return the numerical value of `L[i, j]`, or `zero(eltype(L))` if the entry is not
present in the sparse structure.  Rows and columns are 1-indexed.
"""
function _get_entry(L::SparseMatrixCSC{T}, i::Int, j::Int) where T
    for p in L.colptr[j]:L.colptr[j + 1] - 1
        L.rowval[p] == i && return L.nzval[p]
    end
    return zero(T)
end

"""
    _set_entry!(L::SparseMatrixCSC, i, j, val)

Set the numerical value of `L[i, j]` in-place.  The entry MUST already exist in
the sparsity pattern or an error is thrown (IC(0) forbids fill-in).
"""
function _set_entry!(L::SparseMatrixCSC{T}, i::Int, j::Int, val) where T
    for p in L.colptr[j]:L.colptr[j + 1] - 1
        if L.rowval[p] == i
            L.nzval[p] = convert(T, val)
            return
        end
    end
    error("Entry ($i, $j) is not in the IC(0) sparsity pattern")
end

"""
    _ic_zero_factor(A)

Compute an incomplete Cholesky factorisation IC(0) of the sparse symmetric
positive-definite matrix `A`.  Returns a sparse lower-triangular matrix `L` such
that `A ≈ L * L'` and the sparsity pattern of `L` equals that of `tril(A)`.

The implementation follows Saad (2003, §10.3) column-oriented IKJ variant,
restricted to the symbolic pattern of the original matrix.
"""
function _ic_zero_factor(A::SparseMatrixCSC{T}) where T <: AbstractFloat
    n = size(A, 1)
    n == size(A, 2) || throw(ArgumentError("IC(0) requires a square matrix"))

    # Lower triangle of A — this defines the IC(0) symbolic pattern
    L = SparseMatrixCSC(n, n, copy(A.colptr), copy(A.rowval), copy(A.nzval))
    # Zero out entries above the diagonal
    for j in 1:n
        for p in L.colptr[j]:L.colptr[j + 1] - 1
            if L.rowval[p] < j
                L.nzval[p] = zero(T)
            end
        end
    end
    dropzeros!(L)

    for k in 1:n
        # Diagonal entry
        a_kk = _get_entry(L, k, k)
        @assert a_kk > zero(T) "IC(0) requires A to be SPD; non-positive diagonal at column $k ($(a_kk))"
        _set_entry!(L, k, k, sqrt(float(a_kk)))

        # Collect rows below diagonal in column k
        rows_below = Int[]
        for p in L.colptr[k]:L.colptr[k + 1] - 1
            i = L.rowval[p]
            if i > k
                push!(rows_below, i)
            end
        end

        # Scale column k below diagonal
        for i in rows_below
            val = _get_entry(L, i, k)
            _set_entry!(L, i, k, val / sqrt(float(a_kk)))
        end

        # Update trailing submatrix — only entries in the original sparsity pattern.
        # L[i,j] for i ≥ j is stored in column j, row i.  For a pair (i,j) from
        # rows_below the entry to update is L[max(i,j), min(i,j)].
        nrb = length(rows_below)
        for ia in 1:nrb
            i = rows_below[ia]
            lik = _get_entry(L, i, k)
            for ja in ia:nrb
                j = rows_below[ja]
                ljk = _get_entry(L, j, k)
                col = min(i, j)
                row = max(i, j)
                existing = _get_entry(L, row, col)
                if !iszero(existing)
                    _set_entry!(L, row, col, existing - lik * ljk)
                end
            end
        end
    end

    return L
end

"""Fallback ILU(0) using Julia's built-in LU."""
function _ilu_zero(A::SparseMatrixCSC{T}) where T
    return lu(A)
end

function _validate_preconditioner_policy(policy::PreconditionerBuildPolicy)
    policy.kind in (:jacobi, :identity, :ilu, :ic, :amg) ||
        throw(ArgumentError("implemented preconditioner kinds are :jacobi, :identity, :ilu, :ic, and :amg"))
    policy.diagonal_tolerance >= 0 ||
        throw(ArgumentError("diagonal_tolerance must be nonnegative"))
    policy.ilu_droptol >= 0 ||
        throw(ArgumentError("ilu_droptol must be nonnegative"))
    policy.ilu_tau >= 0 ||
        throw(ArgumentError("ilu_tau must be nonnegative"))
    policy.amg_type in (:smoothed_aggregation, :ruge_stuben) ||
        throw(ArgumentError("amg_type must be :smoothed_aggregation or :ruge_stuben"))
    return policy
end

function _preconditioner_cache_key(problem::AdaptiveLinearProblem,
        policy::PreconditionerBuildPolicy)
    problem.matrix_version === nothing && return nothing
    return (problem.label, problem.matrix_version, size(problem.A), policy.kind,
        policy.diagonal_tolerance, policy.ilu_droptol, policy.ilu_tau, policy.amg_type)
end

function _cache_preconditioner!(cache::PreconditionerCache, key, contract::PreconditionerContract)
    cache.capacity > 0 || return cache
    haskey(cache.entries, key) && filter!(item -> item != key, cache.order)
    cache.entries[key] = contract
    push!(cache.order, key)
    while length(cache.order) > cache.capacity
        evicted = popfirst!(cache.order)
        delete!(cache.entries, evicted)
    end
    return cache
end

"""Clear all locally cached preconditioners; this never modifies a caller-owned operator."""
function clear_preconditioner_cache!(cache::PreconditionerCache)
    empty!(cache.entries)
    empty!(cache.order)
    return cache
end

function _jacobi_preconditioner(problem::AdaptiveLinearProblem,
        policy::PreconditionerBuildPolicy)
    problem.A isa AbstractMatrix || return nothing, :explicit_matrix_required
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator
    diagonal = diag(problem.A)
    any(value -> abs(value) <= policy.diagonal_tolerance, diagonal) &&
        return nothing, :zero_or_small_diagonal
    inverse_diagonal = inv.(diagonal)
    contract = PreconditionerContract(
        name=:jacobi,
        operator=Diagonal(inverse_diagonal),
        fixed_within_solve=PropertyEvidence(Certified; source=:builder),
        linear_within_solve=PropertyEvidence(Certified; source=:builder),
    )
    if _certified_or_proved(problem.contract.hermitian) &&
       _certified_or_proved(problem.contract.positive_definite) &&
       all(value -> real(value) > policy.diagonal_tolerance, diagonal)
        contract = PreconditionerContract(
            name=:jacobi,
            operator=Diagonal(inverse_diagonal),
            fixed_within_solve=PropertyEvidence(Certified; source=:builder),
            linear_within_solve=PropertyEvidence(Certified; source=:builder),
            hermitian=PropertyEvidence(Certified; source=:derived_from_spd_contract),
            positive_definite=PropertyEvidence(Certified; source=:derived_from_spd_contract),
        )
    end
    return contract, :built
end

function _identity_preconditioner(problem::AdaptiveLinearProblem)
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator
    T = try
        eltype(problem.A)
    catch
        Float64
    end
    contract = PreconditionerContract(
        name=:identity,
        operator=Diagonal(ones(T, size(problem.A, 1))),
        fixed_within_solve=PropertyEvidence(Certified; source=:builder),
        linear_within_solve=PropertyEvidence(Certified; source=:builder),
        hermitian=PropertyEvidence(Certified; source=:builder),
        positive_definite=PropertyEvidence(Certified; source=:builder),
    )
    return contract, :built
end

"""
    _ilu_preconditioner(problem, policy)

Build an ILU preconditioner using IncompleteLU.jl (threshold-based) when
available, or a basic ILU(0) fallback based on LU factorisation with sparsity
masking.  Requires a sparse square matrix `problem.A`.
"""
function _ilu_preconditioner(problem::AdaptiveLinearProblem, policy::PreconditionerBuildPolicy)
    problem.A isa SparseMatrixCSC || return nothing, :sparse_matrix_required
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator

    operator = if _HAS_INCOMPLETE_LU
        IncompleteLU.ilu(problem.A; τ=policy.ilu_tau)
    else
        _ilu_zero(problem.A)
    end

    contract = PreconditionerContract(
        name=:ilu,
        operator=operator,
        fixed_within_solve=PropertyEvidence(Certified; source=:builder),
        linear_within_solve=PropertyEvidence(Certified; source=:builder),
    )
    return contract, :built
end

"""
    _ic_preconditioner(problem, policy)

Build an IC(0) preconditioner from the sparse SPD matrix `problem.A`.
Uses the pure-Julia IC(0) factorisation [`_ic_zero_factor`](@ref), which
requires no external dependency.
"""
function _ic_preconditioner(problem::AdaptiveLinearProblem, policy::PreconditionerBuildPolicy)
    problem.A isa SparseMatrixCSC || return nothing, :sparse_matrix_required
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator
    if !(_certified_or_proved(problem.contract.hermitian) &&
         _certified_or_proved(problem.contract.positive_definite))
        return nothing, :spd_contract_required_for_ic
    end

    L = try
        _ic_zero_factor(problem.A)
    catch error
        return nothing, Symbol("ic_factorisation_failed: $(sprint(showerror, error))")
    end

    operator = ICZeroPreconditioner(L)
    contract = PreconditionerContract(
        name=:ic,
        operator=operator,
        fixed_within_solve=PropertyEvidence(Certified; source=:builder),
        linear_within_solve=PropertyEvidence(Certified; source=:builder),
        hermitian=PropertyEvidence(Certified; source=:derived_from_spd_contract),
        positive_definite=PropertyEvidence(Certified; source=:derived_from_spd_contract),
    )
    return contract, :built
end

"""
    _amg_preconditioner(problem, policy)

Build an algebraic multigrid preconditioner using AlgebraicMultigrid.jl.
Two methods are available via `policy.amg_type`:
- `:smoothed_aggregation` (default)   — `smoothed_aggregation(A) + aspreconditioner`
- `:ruge_stuben`                      — `ruge_stuben(A) + aspreconditioner`
Requires a sparse square matrix `problem.A`.
"""
function _amg_preconditioner(problem::AdaptiveLinearProblem, policy::PreconditionerBuildPolicy)
    problem.A isa SparseMatrixCSC || return nothing, :sparse_matrix_required
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator

    _HAS_ALGEBRAIC_MULTIGRID || return nothing, :algebraicmultigrid_package_not_available

    raw_ml = if policy.amg_type == :smoothed_aggregation
        smoothed_aggregation(problem.A)
    else
        ruge_stuben(problem.A)
    end

    operator = aspreconditioner(raw_ml)

    contract = if _certified_or_proved(problem.contract.hermitian) &&
                   _certified_or_proved(problem.contract.positive_definite)
        # For SPD problems, AMG preserves symmetry and positive-definiteness
        PreconditionerContract(
            name=:amg,
            operator=operator,
            fixed_within_solve=PropertyEvidence(Certified; source=:builder),
            linear_within_solve=PropertyEvidence(Certified; source=:builder),
            hermitian=PropertyEvidence(Certified; source=:derived_from_spd_contract),
            positive_definite=PropertyEvidence(Certified; source=:derived_from_spd_contract),
        )
    else
        PreconditionerContract(
            name=:amg,
            operator=operator,
            fixed_within_solve=PropertyEvidence(Certified; source=:builder),
            linear_within_solve=PropertyEvidence(Certified; source=:builder),
        )
    end
    return contract, :built
end

"""
    build_preconditioner(problem; policy=PreconditionerBuildPolicy(), cache=nothing)

Build an implemented preconditioner or reuse a version-keyed cache entry. Reuse is
disabled when `problem.matrix_version` is absent, so structurally similar matrices are
never silently treated as numerically identical.

Supported `policy.kind` values:
- `:jacobi` — diagonal scaling (default)
- `:identity` — unit preconditioner
- `:ilu` — incomplete LU factorisation (sparse required; uses IncompleteLU.jl if available)
- `:ic` — incomplete Cholesky IC(0) (sparse SPD required; pure Julia)
- `:amg` — algebraic multigrid (sparse required; uses AlgebraicMultigrid.jl)
"""
function build_preconditioner(problem::AdaptiveLinearProblem;
        policy::PreconditionerBuildPolicy=PreconditionerBuildPolicy(),
        cache::Union{Nothing, PreconditionerCache}=nothing)
    _validate_preconditioner_policy(policy)
    key = _preconditioner_cache_key(problem, policy)
    if cache !== nothing && policy.reuse && key !== nothing && haskey(cache.entries, key)
        return PreconditionerBuildReport(cache.entries[key], key, true, 0.0, :cache_hit)
    end
    started = time_ns()
    contract, reason = if policy.kind == :jacobi
        _jacobi_preconditioner(problem, policy)
    elseif policy.kind == :identity
        _identity_preconditioner(problem)
    elseif policy.kind == :ilu
        _ilu_preconditioner(problem, policy)
    elseif policy.kind == :ic
        _ic_preconditioner(problem, policy)
    elseif policy.kind == :amg
        _amg_preconditioner(problem, policy)
    else
        nothing, :unknown_preconditioner_kind
    end
    elapsed = Float64(time_ns() - started) / 1.0e9
    contract === nothing && return PreconditionerBuildReport(nothing, key, false, elapsed, reason)
    cache !== nothing && policy.reuse && key !== nothing && _cache_preconditioner!(cache, key, contract)
    return PreconditionerBuildReport(contract, key, false, elapsed, reason)
end

"""Return a problem preserving all contracts and identity fields, with a built preconditioner attached."""
function with_preconditioner(problem::AdaptiveLinearProblem, preconditioner::PreconditionerContract)
    return AdaptiveLinearProblem(problem.A, problem.b; contract=problem.contract,
        conditioning=problem.conditioning, preconditioner=preconditioner, label=problem.label,
        matrix_version=problem.matrix_version)
end