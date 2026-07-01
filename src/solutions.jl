#!/usr/bin/julia

"""
Tools to find voltage solutions (in terms of voltages on the electrodes)
"""
module Solutions

import ..PolyFit, ..gradient
using NLsolve

function find_flat_point(data::A; init=ntuple(i->(size(data, i) + 1) / 2, Val(N))) where (A<:AbstractArray{T,N} where T) where N
    fitter = PolyFit.Fitter(ntuple(i->3, Val(N))...)
    cache = PolyFit.FitCache(fitter, data)
    function model!(g, x)
        xt = ntuple(i->x[i], Val(N))
        for i in 1:N
            g[i] = gradient(cache, i, xt...)
        end
    end
    res = nlsolve(model!, collect(init))
    return ntuple(i->res.zero[i], Val(N))
end

function find_all_flat_points(all_data::A; init=ntuple(i->(size(all_data, i) + 1) / 2, Val(N - 1))) where (A<:AbstractArray{T,N} where T) where N

    npoints = size(all_data, N)
    all_res = Matrix{Float64}(undef, npoints, N - 1)

    for i in 1:npoints
        init = find_flat_point(@view(all_data[:, :, i]), init=init)
        all_res[i, :] .= init
    end
    return all_res
end

struct CenterTracker
    zy_index::Matrix{Float64}
    CenterTracker(zy_index::AbstractMatrix) = new(zy_index)
end

function Base.get(tracker::CenterTracker, xidx)
    # return (y, z)
    nx = size(tracker.zy_index, 1)
    lb_idx = min(max(floor(Int, xidx), 1), nx)
    ub_idx = min(max(ceil(Int, xidx), 1), nx)
    y_lb = tracker.zy_index[lb_idx, 2]
    z_lb = tracker.zy_index[lb_idx, 1]
    if lb_idx == ub_idx
        return y_lb, z_lb
    end
    @assert ub_idx == lb_idx + 1
    y_ub = tracker.zy_index[ub_idx, 2]
    z_ub = tracker.zy_index[ub_idx, 1]
    c_ub = xidx - lb_idx
    c_lb = ub_idx - xidx
    return y_lb * c_lb + y_ub * c_ub, z_lb * c_lb + z_ub * c_ub
end

end
