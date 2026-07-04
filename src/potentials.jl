#!/usr/bin/julia

module Potentials

import ..PolyFit, ..get_single, ..TrapDesc

using DelimitedFiles
using Base.Threads

export Potential, import_pillbox_v0, import_pillbox_v1, import_pillbox_64

public get_potential, Fitting, get_electrodes, load_short_map

struct RawPotential
    electrodes::Int
    nx::Int
    ny::Int
    nz::Int
    stride_um::NTuple{3,Float64}
    origin_um::NTuple{3,Float64}
    data::Array{Float64,4}
end

function import_pillbox_v0_raw(fh)
    read(fh, Int32) # discard
    electrodes = Int(read(fh, Int32))
    nx = Int(read(fh, Int32))
    ny = Int(read(fh, Int32))
    nz = Int(read(fh, Int32))
    read(fh, Int32) # vsets
    stride_um = 1e6 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
    origin_um = 1e6 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
    # I have no idea what's stored in these
    read(fh, Int32)
    read(fh, Int32)
    for i in 1:electrodes
        read(fh, Int32)
    end
    databytes = read(fh)
    if length(databytes) != electrodes * nx * ny * nz * sizeof(Float64)
        throw(ArgumentError("Did not find the right number of samples"))
    end
    data = Array{Float64}(undef, nz, ny, nx, electrodes)
    copyto!(data, reinterpret(Float64, databytes))
    return RawPotential(electrodes, nx, ny, nz, stride_um, origin_um, data)
end

function import_pillbox_v1_raw(fh)
    read(fh, Int32) # discard
    electrodes = Int(read(fh, Int32))
    nx = Int(read(fh, Int32))
    ny = Int(read(fh, Int32))
    nz = Int(read(fh, Int32))
    read(fh, Int32) # vsets
    # xaxis, yaxis
    for _ in 1:6
        read(fh, Float64)
    end
    stride_um = 1e6 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
    origin_um = 1e6 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
    # I have no idea what's stored in these
    read(fh, Int32)
    read(fh, Int32)
    for i in 1:electrodes
        read(fh, Int32)
    end
    databytes = read(fh)
    if length(databytes) != electrodes * nx * ny * nz * sizeof(Float64)
        throw(ArgumentError("Did not find the right number of samples"))
    end
    data = Array{Float64}(undef, nz, ny, nx, electrodes)
    copyto!(data, reinterpret(Float64, databytes))
    return RawPotential(electrodes, nx, ny, nz, stride_um, origin_um, data)
end

function import_pillbox_64_raw(fh)
    read(fh, Int64) # discard
    electrodes = Int(read(fh, Int64))
    nx = Int(read(fh, Int64))
    ny = Int(read(fh, Int64))
    nz = Int(read(fh, Int64))
    read(fh, Int64) # vsets
    # xaxis, yaxis
    for _ in 1:6
        read(fh, Float64)
    end
    stride_um = 1e6 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
    origin_um = 1e6 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
    # I have no idea what's stored in these
    read(fh, Int64)
    read(fh, Int64)
    for i in 1:electrodes
        read(fh, Int64)
    end
    databytes = read(fh)
    if length(databytes) != electrodes * nx * ny * nz * sizeof(Float64)
        throw(ArgumentError("Did not find the right number of samples"))
    end
    data = Array{Float64}(undef, nz, ny, nx, electrodes)
    copyto!(data, reinterpret(Float64, databytes))
    return RawPotential(electrodes, nx, ny, nz, stride_um, origin_um, data)
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        export $(Symbol(name, "_index_to_axis"))
        $(Symbol(name, "_index_to_axis"))(sol::RawPotential, i) = (i - 1) * sol.stride_um[$i] + sol.origin_um[$i]
        export $(Symbol(name, "_axis_to_index"))
        $(Symbol(name, "_axis_to_index"))(sol::RawPotential, a) = (a - sol.origin_um[$i]) / sol.stride_um[$i] + 1
    end
end

const _data_type = typeof(PermutedDimsArray(zeros(0, 0, 0, 0), (3, 2, 1, 4)))

struct Potential
    electrodes::Int
    nx::Int
    ny::Int
    nz::Int
    stride_um::NTuple{3,Float64}
    origin_um::NTuple{3,Float64}
    data::_data_type
    electrode_index::Dict{String,Int}
    electrode_names::Vector{Vector{String}}
    trap::TrapDesc
end

function Potential(raw::RawPotential, electrode_names::AbstractVector,
                   trap::TrapDesc)
    if raw.electrodes != length(trap.ele_names)
        throw(ArgumentError("Electrode number mismatch with trap info."))
    end
    new_electrodes = length(electrode_names)
    data = Array{Float64}(undef, raw.nz, raw.ny, raw.nx, new_electrodes)
    electrode_index = Dict{String,Int}()
    for i in 1:new_electrodes
        electrodes = electrode_names[i]
        first = true
        for elec in electrodes
            electrode_index[elec] = i
            raw_idx = trap.ele_indices[elec]
            if first
                data[:, :, :, i] .= @view raw.data[:, :, :, raw_idx]
                first = false
            else
                data[:, :, :, i] .+= @view raw.data[:, :, :, raw_idx]
            end
        end
        @assert !first
    end
    return Potential(new_electrodes, raw.nx, raw.ny, raw.nz,
                     raw.stride_um, raw.origin_um, PermutedDimsArray(data, (3, 2, 1, 4)),
                     electrode_index, electrode_names, trap)
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        export $(Symbol(name, "_index_to_axis"))
        $(Symbol(name, "_index_to_axis"))(sol::Potential, i) = (i - 1) * sol.stride_um[$i] + sol.origin_um[$i]
        export $(Symbol(name, "_axis_to_index"))
        $(Symbol(name, "_axis_to_index"))(sol::Potential, a) = (a - sol.origin_um[$i]) / sol.stride_um[$i] + 1
    end
end

function _aliases_to_names(aliases::AbstractDict{Int,Int}, trap::TrapDesc)
    # Compute the mapping between id's
    raw_electrodes = length(trap.ele_names)
    id_map = zeros(Int, raw_electrodes)
    id = 0
    new_electrodes = raw_electrodes - length(aliases)
    electrode_names = Vector{Vector{String}}(undef, new_electrodes)
    for i in 1:raw_electrodes
        if i in keys(aliases)
            continue
        end
        id += 1
        id_map[i] = id
        electrode_names[id] = [trap.ele_names[i]]
    end
    @assert new_electrodes == id
    for (k, v) in aliases
        # The user should connect directly to the final one
        if v in keys(aliases)
            throw(ArgumentError("Electrode alias cannot contain alias chain"))
        end
        id = id_map[v]
        @assert id != 0
        name = trap.ele_names[k]
        push!(electrode_names[id], name)
    end
    return electrode_names
end

function _aliases_to_names(aliases::AbstractDict{S1,S2} where {S1<:AbstractString,S2<:AbstractString}, trap::TrapDesc)
    return _aliases_to_names(Dict(trap.ele_indices[k]=>trap.ele_indices[v]
                                  for (k, v) in aliases),
                             trap)
end

function _get_electrode_names(aliases, electrode_names, trap::TrapDesc)
    if electrode_names !== nothing
        if aliases !== nothing
            throw(ArgumentError("At most one of electrode name array or alias map can be probided"))
        end
        return electrode_names
    end
    if aliases === nothing
        return [[name] for name in trap.ele_names]
    end
    return _aliases_to_names(aliases, trap)
end

_read_file(file::AbstractString, @specialize(cb)) = open(cb, file)
_read_file(file, cb) = cb(file)

function _import_internal(file, @specialize(cb), aliases, electrode_names, trap)
    trap = TrapDesc(trap)
    return Potential(_read_file(file, cb),
                     _get_electrode_names(aliases, electrode_names, trap), trap)
end

"""
    import_pillbox_v0(file)::Potential

Imports voltage array files of format V0, here `file` is a `IO` object
or the name of the file to be read.
Voltage arrays contain potential data for one electrode at 1V
and all other electrodes at 0V on a 3D grid of points.

The data is returned in the header and data fields.
header contains the fields 'electrodes' for the number of potentials
for different electrodes,
'nx', 'ny', 'nz' are the number of samples in x-, y-, and z- direction.
'origin_um' is one end point of the 3D grid,
'stride_um' is the stepsize in the 3 dimensions.
"""
function import_pillbox_v0(file; aliases=nothing, electrode_names=nothing, trap)
    return _import_internal(file, import_pillbox_v0_raw, aliases, electrode_names, trap)
end

"""
    import_pillbox_v1(file)::Potential

Imports voltage array files of format V1, the current format,
here `file` is a `IO` object or the name of the file to be read.
Voltage arrays contain potential data for one electrode at 1V
and all other electrodes at 0V on a 3D grid of points.

The data is returned in the header and data fields.
header contains the fields 'electrodes' for the number of potentials
for different electrodes, 'nx', 'ny', 'nz' are the number of samples
in x-, y-, and z- direction.
'origin_um' is one end point of the 3D grid,
'stride_um' is the stepsize in the 3 dimensions.
"""
function import_pillbox_v1(file; aliases=nothing, electrode_names=nothing, trap)
    return _import_internal(file, import_pillbox_v1_raw, aliases, electrode_names, trap)
end

"""
    import_pillbox_64(file)::Potential

Imports voltage array files in 64 bit format,
here `file` is a `IO` object or the name of the file to be read.
Voltage arrays contain potential data for one electrode at 1V
and all other electrodes at 0V on a 3D grid of points.

The data is returned in the header and data fields.
header contains the fields 'electrodes' for the number of potentials
for different electrodes, 'nx', 'ny', 'nz' are the number of samples
in x-, y-, and z- direction.
'origin_um' is one end point of the 3D grid,
'stride_um' is the stepsize in the 3 dimensions.
"""
function import_pillbox_64(file; aliases=nothing, electrode_names=nothing, trap)
    return _import_internal(file, import_pillbox_64_raw, aliases, electrode_names, trap)
end

"""
Calculate the 3D potential given the electrode voltages.
The potential is provided in a `Potential` object whereas the voltages are given
in a map between electrode name/index and the voltage value.
"""
function get_potential(potential::Potential, electrodes_voltages)
    res = zeros(size(potential.data)[1:3])
    for (ele, v) in electrodes_voltages
        if isa(ele, Integer)
            ele_id = Int(ele)::Int
        else
            ele_id = potential.electrode_index[ele]
        end
        res .= muladd.(v, @view(potential.data[:, :, :, ele_id]), res)
    end
    return res
end

const _subarray_T = typeof(@view(PermutedDimsArray(zeros(0, 0, 0, 1), (3, 2, 1, 4))[:, :, :, 1]))

struct Fitting
    fitter::PolyFit.Fitter{3}
    potential::Potential
    cache::AtomicMemory{PolyFit.FitCache{3,_subarray_T}}
    # orders and sizes are in x, y, z order, potential in z, y, x order
    function Fitting(potential::Potential; orders=(4, 2, 2), sizes)
        fitter = PolyFit.Fitter(orders..., sizes=sizes)
        return new(fitter, potential,
                   AtomicMemory{PolyFit.FitCache{3,_subarray_T}}(undef, potential.electrodes))
    end
end

function Base.get(fitting::Fitting, idx::Integer)
    if isassigned(fitting.cache, idx)
        return @atomic :acquire fitting.cache[idx]
    end
    fit_cache = PolyFit.FitCache(fitting.fitter, @view fitting.potential.data[:, :, :, idx])
    cache_ary = fitting.cache
    @atomic :release cache_ary[idx] = fit_cache
    return fit_cache
end

Base.get(fitting::Fitting, name::AbstractString) =
    get(fitting, fitting.potential.electrode_index[name])

Base.get(fitting::Fitting, electrode::Union{AbstractString,Integer},
         pos::NTuple{3}; fit_center=pos) =
             get(get(fitting, electrode), pos; fit_center=fit_center)

get_single(fitting::Fitting, electrode::Union{AbstractString,Integer},
           pos::NTuple{3}, orders::NTuple{3}; fit_center=pos) =
               get_single(get(fitting, electrode), pos, orders;
                          fit_center=fit_center)

function get_electrodes(fitting::Fitting, electrodes_voltages, pos::NTuple{3})
    local res
    for (ele, v) in electrodes_voltages
        term = get(fitting, ele, pos) * v
        if !@isdefined(res)
            res = term
        else
            res += term
        end
    end
    if !@isdefined(res)
        return PolyFit.Result{3}(fitting.fitter.orders,
                                 zeros(prod(fitting.fitter.orders .+ 1)))
    end
    return res
end

function get_electrodes(fitting::Fitting, electrodes_voltages, pos::NTuple{3},
                        orders::NTuple{3})
    res = 0.0
    for (ele, v) in electrodes_voltages
        res += get_single(fitting, ele, pos, orders) * v
    end
    return res
end

function load_short_map(fname)
    m = readdlm(fname, ',', String)
    res = Dict{String,String}()
    for i in 1:size(m, 1)
        res[m[i, 1]] = m[i, 2]
    end
    return res
end

end
