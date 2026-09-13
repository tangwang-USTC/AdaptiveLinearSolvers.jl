"""Mathematical requirements of an iterative method before any numerical kernel is selected."""
struct IterativeMethodCapability
    method::Symbol
    requires_hermitian::Bool
    requires_positive_definite::Bool
    requires_fixed_linear_preconditioner::Bool
    accepts_variable_preconditioner::Bool
end

"""Eligibility result for a requested iterative method; no iteration is performed here."""
struct IterativeQualification
    method::Symbol
    eligible::Bool
    reason::Symbol
    capability::Union{Nothing, IterativeMethodCapability}
end

const _ITERATIVE_METHODS = Dict(
    :cg => IterativeMethodCapability(:cg, true, true, true, false),
    :minres => IterativeMethodCapability(:minres, true, false, true, false),
    :gmres => IterativeMethodCapability(:gmres, false, false, true, false),
    :fgmres => IterativeMethodCapability(:fgmres, false, false, false, true),
    :bicgstab => IterativeMethodCapability(:bicgstab, false, false, true, false),
    :lsqr => IterativeMethodCapability(:lsqr, false, false, true, false),
    :lsmr => IterativeMethodCapability(:lsmr, false, false, true, false),
)

iterative_capability(method::Symbol) = get(_ITERATIVE_METHODS, method, nothing)

function _fixed_linear_preconditioner(problem::AdaptiveLinearProblem)
    preconditioner = problem.preconditioner
    preconditioner === nothing && return true
    return _certified_or_proved(preconditioner.fixed_within_solve) &&
           _certified_or_proved(preconditioner.linear_within_solve)
end

function _hermitian_positive_definite_preconditioner(problem::AdaptiveLinearProblem)
    preconditioner = problem.preconditioner
    preconditioner === nothing && return true
    return _certified_or_proved(preconditioner.hermitian) &&
           _certified_or_proved(preconditioner.positive_definite)
end

"""
    qualify_iterative(problem, method)

Check only mathematical eligibility. In particular, GMRES (Generalized Minimal
Residual) is rejected for an explicitly supplied non-fixed or nonlinear
preconditioner, while FGMRES (Flexible Generalized Minimal Residual) remains eligible.
"""
function qualify_iterative(problem::AdaptiveLinearProblem, method::Symbol)
    capability = iterative_capability(method)
    capability === nothing && return IterativeQualification(method, false, :unknown_iterative_method, nothing)

    contract = problem.contract
    capability.requires_hermitian && !_certified_or_proved(contract.hermitian) &&
        return IterativeQualification(method, false, :hermitian_evidence_insufficient, capability)
    capability.requires_positive_definite && !_certified_or_proved(contract.positive_definite) &&
        return IterativeQualification(method, false, :positive_definite_evidence_insufficient, capability)

    if method in (:cg, :minres) && !_hermitian_positive_definite_preconditioner(problem)
        return IterativeQualification(method, false,
            :preconditioner_hermitian_positive_definite_evidence_insufficient, capability)
    end

    if capability.requires_fixed_linear_preconditioner && !_fixed_linear_preconditioner(problem)
        return IterativeQualification(method, false,
            :variable_or_unknown_preconditioner_requires_fgmres, capability)
    end
    return IterativeQualification(method, true, :qualified, capability)
end
