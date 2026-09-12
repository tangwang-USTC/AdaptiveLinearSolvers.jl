function _iterativesolvers_preconditioner(preconditioner::Union{Nothing, PreconditionerContract})
    (preconditioner === nothing || preconditioner.name == :none) && return NamedTuple()
    preconditioner.operator === nothing &&
        throw(ArgumentError("preconditioner $(preconditioner.name) has no inverse-action operator"))
    return (; Pl=preconditioner.operator)
end

function _iterativesolvers_history(history)
    try
        return Float64[value for value in history[:resnorm]]
    catch
        return Float64[]
    end
end

"""Execute the currently supported IterativeSolvers.jl GMRES adapter."""
function _solve_iterativesolvers(method::Symbol, A, b, residual_policy::ResidualPolicy,
        control::IterationControl, preconditioner::Union{Nothing, PreconditionerContract})
    method == :gmres || throw(ArgumentError("IterativeSolvers adapter only implements GMRES"))
    started = time_ns()
    maximum_iterations = control.max_iterations == 0 ? size(A, 2) : control.max_iterations
    keywords = _iterativesolvers_preconditioner(preconditioner)
    x, history = IterativeSolvers.gmres(A, b;
        reltol=residual_policy.relative_tolerance,
        abstol=residual_policy.absolute_tolerance,
        maxiter=maximum_iterations,
        log=true,
        keywords...)
    elapsed = Float64(time_ns() - started) / 1.0e9
    converged = Bool(getproperty(history, :isconverged))
    iterations = Int(getproperty(history, :iters))
    status = converged ? "converged" : "not converged"
    return x, IterationReport(method, iterations, converged, status, elapsed,
        _iterativesolvers_history(history))
end
