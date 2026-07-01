module Optimizers

import HiGHS
using JuMP
using LinearAlgebra

@inline function gen_minmax_model(B, x0)
    nx = size(B, 1)
    nt = size(B, 2)
    @assert nx == length(x0)

    model = Model(HiGHS.Optimizer)
    set_attribute(model, HiGHS.ComputeInfeasibilityCertificate(), false)
    set_attribute(model, "output_flag", false)
    @variable(model, t[1:nt])
    @variable(model, maxv)
    expr = B * t .+ x0
    @constraint(model, maxv .>= expr)
    @constraint(model, maxv .>= .-expr)
    @objective(model, Min, maxv)
    return model, t
end

@inline function gen_minmax_model_with_limit_terms(B, x0, limited)
    nx = size(B, 1)
    nt = size(B, 2)
    @assert nx == length(x0)
    @assert nx == size(limited, 1)
    nl = size(limited, 2)

    model = Model(HiGHS.Optimizer)
    set_attribute(model, HiGHS.ComputeInfeasibilityCertificate(), false)
    set_attribute(model, "output_flag", false)
    @variable(model, t[1:nt])
    @variable(model, l[1:nl])
    @variable(model, maxv)
    expr = B * t .+ x0 .+ limited * l
    @constraint(model, maxv .>= expr)
    @constraint(model, maxv .>= .-expr)
    @constraint(model, l .<= 1)
    @constraint(model, l .>= -1)
    @objective(model, Min, maxv)
    return model, t, l
end

function optimize_minmax_span(B, x0; limited=nothing)
    if limited === nothing
        model, t = gen_minmax_model(B, x0)
        JuMP.optimize!(model)
        return B * value.(t) .+ x0
    else
        model, t, l = gen_minmax_model_with_limit_terms(B, x0, limited)
        JuMP.optimize!(model)
        return B * value.(t) .+ x0 .+ limited * value.(l)
    end
end

function optimize_minmax(A, y::Union{AbstractVector, AbstractMatrix})
    x0 = A \ y
    ny, nx = size(A)
    ns = size(y, 2)
    nt = nx - ny
    if nt <= 0
        return x0
    end
    # With the A * x = y constraints,
    # the degrees of freedom left in x are the ones that satisfies A * x = 0
    # In another word, these are the x's that are orthogonal to all rows of A.
    # We can find the basis set that spans such space using QR decomposition.
    B = @view(qr(A').Q[:, (ny + 1):nx])
    for i in 1:ns
        x = @view(x0[:, i])
        model, t = gen_minmax_model(B, x)
        JuMP.optimize!(model)
        mul!(x, B, value.(t), true, true)
    end
    return x0
end

end
