function _validate_preconditioner_policy(policy::PreconditionerBuildPolicy)
    policy.kind in (:jacobi, :identity) ||
        throw(ArgumentError("implemented preconditioner kinds are :jacobi and :identity"))
    policy.diagonal_tolerance >= 0 ||
        throw(ArgumentError("diagonal_tolerance must be nonnegative"))
    return policy
end

function _preconditioner_cache_key(problem::AdaptiveLinearProblem,
        policy::PreconditionerBuildPolicy)
    problem.matrix_version === nothing && return nothing
    return (problem.label, problem.matrix_version, size(problem.A), policy.kind,
        policy.diagonal_tolerance)
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
    build_preconditioner(problem; policy=PreconditionerBuildPolicy(), cache=nothing)

Build an implemented preconditioner or reuse a version-keyed cache entry. Reuse is
disabled when `problem.matrix_version` is absent, so structurally similar matrices are
never silently treated as numerically identical.
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
    contract, reason = policy.kind == :jacobi ? _jacobi_preconditioner(problem, policy) :
        _identity_preconditioner(problem)
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
