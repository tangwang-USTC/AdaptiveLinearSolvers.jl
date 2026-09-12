const _CONDITIONING_METRICS = (:kappa_1, :kappa_2, :kappa_inf, :rcond_1, :rcond_2, :rcond_inf)
const _CONDITIONING_OPERATORS = (:original, :left_preconditioned, :right_preconditioned)
const _CONDITIONING_EVIDENCE = (:exact, :estimated, :external, :qualitative)

function _validate_conditioning_policy(policy::ConditioningPolicy)
    policy.estimation in (:none, :cheap, :full) ||
        throw(ArgumentError("conditioning estimation must be :none, :cheap, or :full"))
    policy.spectral_estimation in (:none, :lanczos) ||
        throw(ArgumentError("spectral_estimation must be :none or :lanczos"))
    policy.spectral_steps > 0 || throw(ArgumentError("spectral_steps must be positive"))
    policy.budget.max_seconds >= 0 ||
        throw(ArgumentError("diagnostic max_seconds must be nonnegative"))
    policy.budget.max_operator_applications >= 0 ||
        throw(ArgumentError("diagnostic max_operator_applications must be nonnegative"))
    policy.budget.max_matrix_dimension >= 0 ||
        throw(ArgumentError("diagnostic max_matrix_dimension must be nonnegative"))
    return policy
end

function _lanczos_spectrum(problem::AdaptiveLinearProblem, policy::ConditioningPolicy)
    policy.spectral_estimation == :none && return nothing, :spectral_estimation_disabled, 0.0
    _certified_or_proved(problem.contract.hermitian) ||
        return nothing, :hermitian_evidence_insufficient, 0.0
    budget = policy.budget
    budget.max_seconds > 0 || return nothing, :diagnostic_time_budget_not_granted, 0.0
    budget.max_operator_applications > 0 ||
        return nothing, :operator_application_budget_not_granted, 0.0
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator, 0.0
    steps = min(policy.spectral_steps, budget.max_operator_applications, size(problem.A, 1))
    steps > 0 || return nothing, :empty_spectral_budget, 0.0

    T = try
        promote_type(eltype(problem.A), Float64)
    catch
        Float64
    end
    q_previous = zeros(T, size(problem.A, 2))
    q = ones(T, size(problem.A, 2))
    q ./= norm(q)
    work = similar(q)
    diagonal = Float64[]
    offdiagonal = Float64[]
    beta_previous = 0.0
    started = time_ns()
    for step in 1:steps
        mul!(work, problem.A, q)
        alpha = Float64(real(dot(q, work)))
        work .-= alpha .* q
        step > 1 && (work .-= beta_previous .* q_previous)
        beta = Float64(norm(work))
        push!(diagonal, alpha)
        if step == steps || beta <= sqrt(eps(Float64))
            break
        end
        push!(offdiagonal, beta)
        q_previous .= q
        q .= work ./ beta
        beta_previous = beta
        elapsed = Float64(time_ns() - started) / 1.0e9
        elapsed <= budget.max_seconds ||
            return nothing, :diagnostic_time_budget_exceeded, elapsed
    end
    elapsed = Float64(time_ns() - started) / 1.0e9
    elapsed <= budget.max_seconds || return nothing, :diagnostic_time_budget_exceeded, elapsed
    values = eigvals(SymTridiagonal(diagonal, offdiagonal))
    info = SpectralInfo(lambda_min=Float64(minimum(values)), lambda_max=Float64(maximum(values)),
        method=:lanczos, operator=:original, matrix_version=problem.matrix_version,
        source=:router_estimate, steps=length(diagonal), reliable=false)
    return info, :estimated, elapsed
end

function _spectral_diagnostic_state(problem::AdaptiveLinearProblem,
        information::Union{Nothing, SpectralInfo})
    information === nothing && return :unavailable
    _certified_or_proved(problem.contract.positive_definite) || return :estimated
    scale = max(abs(information.lambda_max), 1.0)
    information.lambda_min <= sqrt(eps(Float64)) * scale && return :near_singular
    return :estimated
end

function _condition_number(info::ConditioningInfo)
    info.estimate === nothing && return nothing
    return info.metric in (:rcond_1, :rcond_2, :rcond_inf) ?
        (info.estimate == 0 ? Inf : inv(info.estimate)) : info.estimate
end

function _accepted_conditioning_state(problem::AdaptiveLinearProblem, info::ConditioningInfo)
    _certified_or_proved(problem.contract.rank_deficient) && return :near_rank_deficient
    kappa = _condition_number(info)
    kappa === nothing && return :unknown
    kappa >= inv(sqrt(eps(Float64))) && return :near_rank_deficient
    kappa >= 1e6 && return :ill_conditioned
    return :moderate
end

function _assess_conditioning(problem::AdaptiveLinearProblem, info::Union{Nothing, ConditioningInfo},
        policy::ConditioningPolicy; generated::Bool=false)
    info === nothing && return ConditioningAssessment(nothing, false, :no_conditioning_information, :unknown)
    !generated && !policy.use_external &&
        return ConditioningAssessment(info, false, :external_conditioning_disabled, :unknown)
    info.metric in _CONDITIONING_METRICS ||
        return ConditioningAssessment(info, false, :unsupported_conditioning_metric, :invalid)
    info.operator in _CONDITIONING_OPERATORS ||
        return ConditioningAssessment(info, false, :unsupported_conditioning_operator, :invalid)
    info.evidence in _CONDITIONING_EVIDENCE ||
        return ConditioningAssessment(info, false, :unsupported_conditioning_evidence, :invalid)
    (info.estimate === nothing || info.estimate < 0 || isnan(info.estimate)) &&
        return ConditioningAssessment(info, false, :invalid_conditioning_estimate, :invalid)
    policy.require_matching_version &&
        (problem.matrix_version === nothing || info.matrix_version === nothing ||
         info.matrix_version != problem.matrix_version) &&
        return ConditioningAssessment(info, false, :matrix_version_mismatch, :stale)
    !generated && policy.require_reliable_external && !info.reliable &&
        return ConditioningAssessment(info, false, :external_conditioning_not_reliable, :untrusted)
    return ConditioningAssessment(info, true, :accepted, _accepted_conditioning_state(problem, info))
end

"""
    assess_conditioning(problem, policy=ConditioningPolicy())

Validate caller-provided conditioning data without factorizing or applying the operator.
"""
function assess_conditioning(problem::AdaptiveLinearProblem,
        policy::ConditioningPolicy=ConditioningPolicy())
    _validate_conditioning_policy(policy)
    return _assess_conditioning(problem, problem.conditioning, policy)
end

function _estimate_conditioning(problem::AdaptiveLinearProblem, policy::ConditioningPolicy)
    policy.estimation == :none && return nothing, :estimation_disabled, 0.0
    budget = policy.budget
    budget.max_seconds > 0 || return nothing, :diagnostic_time_budget_not_granted, 0.0
    problem.A isa StridedMatrix || return nothing, :explicit_dense_matrix_required, 0.0
    size(problem.A, 1) == size(problem.A, 2) || return nothing, :nonsquare_operator, 0.0
    size(problem.A, 1) <= budget.max_matrix_dimension ||
        return nothing, :diagnostic_matrix_dimension_budget_exceeded, 0.0

    metric = policy.estimation == :cheap ? :rcond_1 : :kappa_2
    started = time_ns()
    value = try
        policy.estimation == :cheap ? inv(cond(problem.A, 1)) : cond(problem.A, 2)
    catch
        return nothing, :conditioning_estimation_failed, Float64(time_ns() - started) / 1.0e9
    end
    elapsed = Float64(time_ns() - started) / 1.0e9
    elapsed <= budget.max_seconds || return nothing, :diagnostic_time_budget_exceeded, elapsed
    isfinite(value) || return nothing, :nonfinite_conditioning_estimate, elapsed
    info = ConditioningInfo(estimate=Float64(value), metric=metric, operator=:original,
        matrix_version=problem.matrix_version, source=:router_estimate, reliable=true,
        evidence=:estimated)
    return info, :estimated, elapsed
end

function _iteration_diagnostic_state(iteration)
    iteration === nothing && return :not_observed
    iteration.converged && return :converged
    backend_status = lowercase(iteration.backend_status)
    occursin("breakdown", backend_status) && return :breakdown
    history = iteration.residual_history
    length(history) >= 3 && history[1] > 0 && history[end] / history[1] >= 0.9 &&
        return :stagnated
    return :nonconverged
end

function _preconditioner_diagnostic_state(problem::AdaptiveLinearProblem, iteration_state::Symbol)
    preconditioner = problem.preconditioner
    (preconditioner === nothing || preconditioner.name == :none) && return :not_used
    iteration_state in (:stagnated, :breakdown) && return :suspected_failure
    return :not_assessed
end

function _overall_diagnostic_state(conditioning::ConditioningAssessment,
        spectral_state::Symbol, iteration_state::Symbol, preconditioner_state::Symbol, status)
    conditioning.state == :near_rank_deficient && return :near_rank_deficient
    spectral_state == :near_singular && return :near_rank_deficient
    iteration_state == :breakdown && return :iterative_breakdown
    iteration_state == :stagnated && return :iterative_stagnation
    preconditioner_state == :suspected_failure && return :preconditioner_suspected_failure
    (iteration_state == :nonconverged || status == NumericalFailure) && return :nonconverged
    conditioning.state == :ill_conditioned && return :ill_conditioned
    return :inconclusive
end

"""
    diagnose(problem; policy=ConditioningPolicy(), iteration=nothing, status=nothing)

Return a budgeted numerical diagnosis. The default policy only validates supplied data;
it never computes a condition estimate or performs an additional operator application.
"""
function diagnose(problem::AdaptiveLinearProblem; policy::ConditioningPolicy=ConditioningPolicy(),
        iteration=nothing, status=nothing)
    _validate_conditioning_policy(policy)
    assessment = assess_conditioning(problem, policy)
    estimate_performed = false
    elapsed = 0.0
    notes = String[]
    if !assessment.accepted && policy.estimation != :none
        information, reason, elapsed = _estimate_conditioning(problem, policy)
        if information === nothing
            push!(notes, string(reason))
        else
            assessment = _assess_conditioning(problem, information, policy; generated=true)
            estimate_performed = assessment.accepted
            push!(notes, string(reason))
        end
    end
    spectral, spectral_reason, spectral_elapsed = _lanczos_spectrum(problem, policy)
    spectral_state = policy.spectral_estimation == :none ? :not_requested :
        _spectral_diagnostic_state(problem, spectral)
    elapsed += spectral_elapsed
    policy.spectral_estimation != :none && push!(notes, string(spectral_reason))
    iteration_state = _iteration_diagnostic_state(iteration)
    preconditioner_state = _preconditioner_diagnostic_state(problem, iteration_state)
    overall_state = _overall_diagnostic_state(assessment, spectral_state, iteration_state,
        preconditioner_state, status)
    return NumericalDiagnosis(assessment, spectral, spectral_state, iteration_state, preconditioner_state,
        overall_state, estimate_performed, elapsed, notes)
end
