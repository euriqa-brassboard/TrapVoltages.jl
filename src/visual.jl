#

module Visual

import ..load_optdep, ..TrapDesc

export SVGInfo
public set_title!, add_circle!, add_line!, fill_electrodes!

mutable struct Title
    title::String
    const attrs::Dict{String,String}
    Title() = new("", Dict{String,String}())
end

Base.getindex(t::Title, name::AbstractString) = t.attrs[name]
function Base.setindex!(t::Title, v, name::AbstractString)
    t.attrs[name] = v
    return
end

struct Circle
    x::Float64
    y::Float64
    r::Float64
    attrs::Dict{String,String}
end

Base.getindex(c::Circle, name::AbstractString) = c.attrs[name]
function Base.setindex!(c::Circle, v, name::AbstractString)
    c.attrs[name] = v
    return
end

struct Line
    xs::Vector{Float64}
    ys::Vector{Float64}
    attrs::Dict{String,String}
end

Base.getindex(l::Line, name::AbstractString) = l.attrs[name]
function Base.setindex!(l::Line, v, name::AbstractString)
    l.attrs[name] = v
    return
end

mutable struct SVGInfo
    const trap::TrapDesc
    const title::Title
    const circles::Vector{Circle}
    const lines::Vector{Line}
    const electrode_colors::Dict{String,String}
    function SVGInfo(trap)
        return new(TrapDesc(trap), Title(), Circle[], Line[], Dict{String,String}())
    end
end

function set_title!(info::SVGInfo, title)
    t = info.title
    t.title = title
    return t
end

function add_circle!(info::SVGInfo, x, y, r)
    c = Circle(x, y, r, Dict{String,String}())
    push!(info.circles, c)
    return c
end

function add_line!(info::SVGInfo, xs, ys)
    l = Line(xs, ys, Dict{String,String}())
    push!(info.lines, l)
    return l
end

const rdbu_cmap_data = ((5, 48, 97), (33, 102, 172), (67, 147, 195),
                        (146, 197, 222), (209, 229, 240),
                        (247, 247, 247),
                        (253, 219, 199), (244, 165, 130), (214, 96, 77),
                        (178, 24, 43), (103, 0, 31))

function default_cmap(v) # v ∈ [-1, 1]
    npts = length(rdbu_cmap_data)
    nsegs = npts - 1
    x = nsegs * (Float64(v) + 1) * 0.5 + 1
    if !(x > 1)
        return rdbu_cmap_data[1]
    elseif !(x < npts)
        return rdbu_cmap_data[npts]
    end
    idx = unsafe_trunc(Int, x)
    x = x - idx
    c1 = rdbu_cmap_data[idx]
    c2 = rdbu_cmap_data[idx + 1]
    return unsafe_trunc.(Int, c2 .* x .+ c1 .* (1 - x) .+ 0.5)
end

const hex_str_data = ((string(i, base=16, pad=2) for i in 0:255)...,)

@inline function hex_str(v::Integer)
    if v < 0
        return @inbounds hex_str_data[1]
    elseif v > 255
        return @inbounds hex_str_data[256]
    else
        return @inbounds hex_str_data[Int(v + 1)]
    end
end

function rgb_str(v::NTuple{3,Integer})
    r, g, b = hex_str.(v)
    return "#" * r * g * b
end

rgb_str(v::NTuple{3,AbstractFloat}) = rgb_str(round.(Int, v .* 255))
rgb_str(v::AbstractString) = String(v)

function _render_ezxml end
const _ezxml_init = Ref(false)
const _ezxml_pkgid = Base.PkgId(Base.UUID("8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"),
                                "EzXML")

function Base.show(io::IO, ::MIME"image/svg+xml", info::SVGInfo)
    load_optdep(_ezxml_pkgid, _ezxml_init, "rendering SVG")
    invokelatest(_render_ezxml, io, info)
    return
end

function fill_electrodes!(info::SVGInfo, values; cmap=default_cmap)
    for (ele, val) in values
        if isa(val, Number)
            val = cmap(val)
        end
        val = rgb_str(val)::String
        info.electrode_colors[ele] = val
    end
end

end
