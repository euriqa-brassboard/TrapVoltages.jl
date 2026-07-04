#!/usr/bin/julia

"""
Tools to find voltage solutions (in terms of voltages on the electrodes)
"""
module Solutions

import ..PolyFit, ..gradient, ..Units.TrapUnits, ..Potentials, ..Optimizers
using ..Traps, ..Units

using NLsolve
using LinearAlgebra

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

    npoints = size(all_data, 1)
    all_res = Matrix{Float64}(undef, npoints, N - 1)

    for i in 1:npoints
        init = find_flat_point(@view(all_data[i, :, :]), init=init)
        all_res[i, :] .= init
    end
    return all_res
end

struct CenterTracker
    yz_index::Matrix{Float64}
    CenterTracker(yz_index::AbstractMatrix) = new(yz_index)
end

function _get_rf_center_path(name)
    if !('/' in name)
        p = joinpath(@__DIR__, "../data", name, "rf_center.h5")
        if isfile(p)
            return p
        end
    end
    return name
end

function _center_tracker_hdf5 end

const _hdf5_init = Ref(false)

function CenterTracker(name::AbstractString, region=1)
    path = _get_rf_center_path(name)
    if !_hdf5_init[]
        Base.require(Base.PkgId(Base.UUID("f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"), "HDF5"))
    end
    _hdf5_init[] = true
    return invokelatest(_center_tracker_hdf5, path, "$region")::CenterTracker
end

function Base.get(tracker::CenterTracker, xidx)
    # return (y, z)
    nx = size(tracker.yz_index, 1)
    lb_idx = min(max(floor(Int, xidx), 1), nx)
    ub_idx = min(max(ceil(Int, xidx), 1), nx)
    y_lb = tracker.yz_index[lb_idx, 1]
    z_lb = tracker.yz_index[lb_idx, 2]
    if lb_idx == ub_idx
        return y_lb, z_lb
    end
    @assert ub_idx == lb_idx + 1
    y_ub = tracker.yz_index[ub_idx, 1]
    z_ub = tracker.yz_index[ub_idx, 2]
    c_ub = xidx - lb_idx
    c_lb = ub_idx - xidx
    return y_lb * c_lb + y_ub * c_ub, z_lb * c_lb + z_ub * c_ub
end

struct TermMask
    dx::Bool # x
    dy::Bool # y
    dz::Bool # z

    xy::Bool # xy
    yz::Bool # yz
    zx::Bool # zx
    z2::Bool # (z^2 - y^2) / 2

    x2::Bool # (x^2 - (y^2 + z^2) / 2) / 2
    x3::Bool # x^3 / 3!
    x4::Bool # x^4 / 4!
    x2z::Bool # x^2z / 2
    global @inline _term_mask(args...) = new(args...)
end

Base.count(mask::TermMask) = (mask.dx + mask.dy + mask.dz + mask.xy + mask.yz + mask.zx +
    mask.z2 + mask.x2 + mask.x3 + mask.x4 + mask.x2z)

@inline TermMask(; dx=true, dy=true, dz=true, xy=true, yz=true, zx=true, z2=true,
                 x2=true, x3=true, x4=true, x2z=false) =
                     Val(_term_mask(dx, dy, dz, xy, yz, zx, z2, x2, x3, x4, x2z))

function compensate_terms(fit::PolyFit.Result{3}, stride_um;
                          unit::TrapUnits, mask::Val{Terms}=TermMask()) where Terms
    # axis order of fitting result is z, y, x
    # axis order of stride_um is x, y, z
    res = (;)

    # Expected units
    # DX/DY/DZ: V/m
    # Rest: follow TrapUnits
    scale_1 = 1e6
    scale_2 = (unit.l_unit_um^2 / unit.V_unit)
    scale_3 = (unit.l_unit_um^3 / unit.V_unit)
    scale_4 = (unit.l_unit_um^4 / unit.V_unit)

    if Terms.dx
        res = (; res..., dx=fit[1, 0, 0] / stride_um[1] * scale_1)
    end
    if Terms.dy
        res = (; res..., dy=fit[0, 1, 0] / stride_um[2] * scale_1)
    end
    if Terms.dz
        res = (; res..., dz=fit[0, 0, 1] / stride_um[3] * scale_1)
    end

    if Terms.xy
        res = (; res..., xy=fit[1, 1, 0] / stride_um[1] / stride_um[2] * scale_2)
    end
    if Terms.yz
        res = (; res..., yz=fit[0, 1, 1] / stride_um[2] / stride_um[3] * scale_2)
    end
    if Terms.zx
        res = (; res..., zx=fit[1, 0, 1] / stride_um[3] / stride_um[1] * scale_2)
    end

    # The two legal quadratic terms are `x^2 - (y^2 + z^2) / 2` and `z^2 - y^2`
    # which are also orthogonal to each other.
    # The orthogonal illegal term is `x^2 + y^2 + z^2`.
    # Here we just need to find the transfermation to go from the taylor expansion
    # basis to the new basis.
    # Since the three terms are orthogonal, we can just compute the dot product
    # with these three terms and apply the correct normalization coefficient.

    # We need to divide the xy/yz/zx terms by 2 relative to the x2, y2, z2 terms.
    # This makes sure that, e.g., the z^2 term is a direct rotation of the
    # xy/yz/zx terms.
    if Terms.x2 || Terms.z2
        y2 = fit[0, 2, 0] / stride_um[2]^2 * 2
        z2 = fit[0, 0, 2] / stride_um[3]^2 * 2
    end

    if Terms.z2
        res = (; res..., z2=(z2 - y2) / 2 * scale_2)
    end
    if Terms.x2
        x2 = fit[2, 0, 0] / stride_um[1]^2 * 2
        res = (; res..., x2=(2 * x2 - y2 - z2) / 3 * scale_2)
    end
    if Terms.x3
        res = (; res..., x3=fit[3, 0, 0] / stride_um[1]^3 * 6 * scale_3)
    end
    if Terms.x4
        res = (; res..., x4=fit[4, 0, 0] / stride_um[1]^4 * 24 * scale_4)
    end
    if Terms.x2z
        res = (; res..., x2z=fit[2, 0, 1] / stride_um[1]^2 / stride_um[3] * 2 * scale_3)
    end

    return res
end

function compensate_terms(fitting::Potentials.Fitting, pos::NTuple{3};
                          unit::TrapUnits, mask::Val{Terms}=TermMask(),
                          min_num=20, min_dist=0, region=1, ignore_id=(1,)) where Terms
    # pos is in xyz index
    potential = fitting.potential
    x_um = Potentials.x_index_to_axis(potential, pos[1])
    ele_select = find_electrodes(potential.trap, potential.electrode_index,
                                 x_um, min_num=min_num, min_dist=min_dist,
                                 region=1, ignore_id=ignore_id)
    ele_select = sort!(collect(ele_select))
    fits = [get(fitting, e, pos) for e in ele_select]
    stride_um = potential.stride_um
    nfits = length(fits)
    nterms = count(Terms)
    coefficient = Matrix{Float64}(undef, nterms, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(compensate_terms(fits[i], stride_um,
                                                    unit=unit, mask=mask))
    end
    return ele_select, coefficient
end

function solve_compensate(fitting::Potentials.Fitting, pos::NTuple{3};
                          mask::Val{Terms}=TermMask(), minmax=true, kws...) where Terms
    ele_select, coefficient = compensate_terms(fitting, pos; mask=mask, kws...)
    nterms = count(Terms)
    M = Matrix(I, nterms, nterms)
    if minmax
        X = Optimizers.optimize_minmax(coefficient, M)
    else
        X = coefficient \ M
    end
    res = (;)
    idx = 1
    if Terms.dx
        res = (; res..., dx=X[:, idx])
        idx += 1
    end
    if Terms.dy
        res = (; res..., dy=X[:, idx])
        idx += 1
    end
    if Terms.dz
        res = (; res..., dz=X[:, idx])
        idx += 1
    end
    if Terms.xy
        res = (; res..., xy=X[:, idx])
        idx += 1
    end
    if Terms.yz
        res = (; res..., yz=X[:, idx])
        idx += 1
    end
    if Terms.zx
        res = (; res..., zx=X[:, idx])
        idx += 1
    end
    if Terms.z2
        res = (; res..., z2=X[:, idx])
        idx += 1
    end
    if Terms.x2
        res = (; res..., x2=X[:, idx])
        idx += 1
    end
    if Terms.x3
        res = (; res..., x3=X[:, idx])
        idx += 1
    end
    if Terms.x4
        res = (; res..., x4=X[:, idx])
        idx += 1
    end
    if Terms.x2z
        res = (; res..., x2z=X[:, idx])
        idx += 1
    end
    return ele_select, res
end

function solve_target(fitting::Potentials.Fitting, pos::NTuple{3}, target;
                      mask::Val{Terms}=TermMask(), minmax=true, kws...) where Terms
    nterms = count(Terms)
    if size(target, 1) != nterms
        throw(ArgumentError("Target size mismatch with term number"))
    end
    ele_select, coefficient = compensate_terms(fitting, pos; mask=mask, kws...)
    if minmax
        return ele_select, Optimizers.optimize_minmax(coefficient, target)
    else
        return ele_select, coefficient \ target
    end
end

function __init__()
    _hdf5_init[] = false
end

end
