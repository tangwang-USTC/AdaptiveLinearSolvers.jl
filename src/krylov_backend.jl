"""Normalized result from one Krylov.jl invocation."""
struct IterationReport
    method::Symbol
    iterations::Int
    converged::Bool
    backend_status::String
    elapsed_seconds::Float64
    residual_history::Vector{Float64}
    operator_applications::Union{Nothing, Int}
end

IterationReport(method::Symbol, iterations::Int, converged::Bool, backend_status::String,
    elapsed_seconds::Float64, residual_history::Vector{Float64}) =
    IterationReport(method, iterations, converged, backend_status, elapsed_seconds,
        residual_history, nothing)

function _validate_iteration_control(control::IterationControl)
    control.max_iterations >= 0 || throw(ArgumentError("max_iterations must be nonnegative"))
    control.max_seconds > 0 || throw(ArgumentError("max_seconds must be positive"))
    return control
end

function _krylov_history(stats)
    hasproperty(stats, :residuals) || return Float64[]
    return Float64[value for value in getproperty(stats, :residuals)]
end

function _krylov_report(method::Symbol, stats, started_ns::UInt64)
    elapsed_seconds = Float64(time_ns() - started_ns) / 1.0e9
    return IterationReport(method, Int(getproperty(stats, :niter)),
        Bool(getproperty(stats, :solved)), String(getproperty(stats, :status)),
        elapsed_seconds, _krylov_history(stats))
end

function _solve_krylov(method::Symbol, A, b, residual_policy::ResidualPolicy,
        control::IterationControl, preconditioner::Union{Nothing, PreconditionerContract})
    started_ns = time_ns()
    common = (; atol=residual_policy.absolute_tolerance,
        rtol=residual_policy.relative_tolerance,
        itmax=control.max_iterations,
        timemax=control.max_seconds,
        verbose=0,
        history=control.record_history)
    preconditioner_keywords = _krylov_preconditioner_keywords(method, preconditioner)
    x, stats = if method == :cg
        Krylov.cg(A, b; common..., preconditioner_keywords...)
    elseif method == :minres
        Krylov.minres(A, b; common..., preconditioner_keywords...)
    elseif method == :gmres
        Krylov.gmres(A, b; common..., preconditioner_keywords..., restart=control.restart)
    elseif method == :fgmres
        Krylov.fgmres(A, b; common..., preconditioner_keywords..., restart=control.restart)
    elseif method == :bicgstab
        Krylov.bicgstab(A, b; common..., preconditioner_keywords...)
    else
        throw(ArgumentError("Krylov backend has no implementation for route $method"))
    end
    return x, _krylov_report(method, stats, started_ns)
end

function _krylov_preconditioner_keywords(method::Symbol,
        preconditioner::Union{Nothing, PreconditionerContract})
    (preconditioner === nothing || preconditioner.name == :none) && return NamedTuple()
    preconditioner.operator === nothing &&
        throw(ArgumentError("preconditioner $(preconditioner.name) has no inverse-action operator"))
    return method == :fgmres ? (; N=preconditioner.operator) : (; M=preconditioner.operator)
end

function _iterative_backend_available(problem::AdaptiveLinearProblem)
    preconditioner = problem.preconditioner
    return preconditioner === nothing || preconditioner.name == :none ||
           preconditioner.operator !== nothing
end

function _iteration_failure_status(report::IterationReport, control::IterationControl)
    status = lowercase(report.backend_status)
    if (control.max_iterations > 0 && report.iterations >= control.max_iterations) ||
       (isfinite(control.max_seconds) && report.elapsed_seconds >= control.max_seconds) ||
       occursin("time limit", status) || occursin("maximum iterations", status)
        return BudgetTerminated
    end
    return NumericalFailure
end
