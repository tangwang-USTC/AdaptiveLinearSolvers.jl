function _linearsolve_success(solution)
    try
        return Bool(LinearSolve.SciMLBase.successful_retcode(solution.retcode))
    catch
        return false
    end
end

"""Execute the deliberately narrow LinearSolve.jl GMRES adapter."""
function _solve_linearsolve(method::Symbol, A, b, residual_policy::ResidualPolicy,
        control::IterationControl, preconditioner::Union{Nothing, PreconditionerContract})
    method == :gmres || throw(ArgumentError("LinearSolve adapter only implements GMRES"))
    A isa AbstractMatrix || throw(ArgumentError("LinearSolve adapter currently requires an explicit matrix"))
    (preconditioner === nothing || preconditioner.name == :none) ||
        throw(ArgumentError("LinearSolve adapter does not yet map user preconditioners"))
    started = time_ns()
    problem = LinearSolve.LinearProblem(A, b)
    solution = LinearSolve.solve(problem, LinearSolve.KrylovJL_GMRES();
        abstol=residual_policy.absolute_tolerance,
        reltol=residual_policy.relative_tolerance,
        maxiters=control.max_iterations == 0 ? size(A, 2) : control.max_iterations)
    elapsed = Float64(time_ns() - started) / 1.0e9
    converged = _linearsolve_success(solution)
    status = string(solution.retcode)
    return solution.u, IterationReport(method, 0, converged, status, elapsed, Float64[])
end
