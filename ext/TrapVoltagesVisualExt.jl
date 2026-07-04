#

module TrapVoltagesVisualExt

using TrapVoltages.Visual
using EzXML
const svg_ns = ["svg" => "http://www.w3.org/2000/svg"]

function get_template(trap_name, has_title, has_plot)
    if has_plot
        name = "extra-space.svg"
    elseif has_title
        name = "title.svg"
    else
        name = "raw.svg"
    end
    return readxml(joinpath(@__DIR__, "../data", trap_name, "templates", name))
end

function apply_attrs!(node, attrs)
    for (k, v) in attrs
        node[k] = v
    end
end

struct SVGCoord
    # The transformation is
    # sx = tx * C + ty * S + X
    # sy = -tx * S + ty * C + Y
    C::Float64
    S::Float64
    X::Float64
    Y::Float64
    scale::Float64
end

function SVGCoord(svg)
    ionpos = findfirst("//svg:line[@id='IonPos']", svg, svg_ns)
    @assert ionpos !== nothing
    sx1 = parse(Float64, ionpos["x1"])
    sy1 = parse(Float64, ionpos["y1"])
    sx2 = parse(Float64, ionpos["x2"])
    sy2 = parse(Float64, ionpos["y2"])

    tx1 = parse(Float64, ionpos["data-trap-x1"])
    ty1 = parse(Float64, ionpos["data-trap-y1"])
    tx2 = parse(Float64, ionpos["data-trap-x2"])
    ty2 = parse(Float64, ionpos["data-trap-y2"])
    unlink!(ionpos)

    (C, S, X, Y) = [tx1 ty1 1 0
                    ty1 -tx1 0 1
                    tx2 ty2 1 0
                    ty2 -tx2 0 1] \ [sx1, sy1, sx2, sy2]
    return SVGCoord(C, S, X, Y, hypot(C, S))
end

to_svg(sc::SVGCoord, tx, ty) =
    (tx * sc.C + ty * sc.S + sc.X, -tx * sc.S + ty * sc.C + sc.Y)
to_svg(sc::SVGCoord, r) = sc.scale * r

function num2str(v)
    v = round(v, digits=3)
    if isinteger(v)
        return string(Int(v))
    end
    return string(v)
end

function draw_line(svg, l, coord, ytop, height)
    # This currently assume the plotting area is horizontal
    line = addelement!(svg, "polyline")
    points = String[]
    local last_y, pending
    is_pending = false
    for (x, y) in zip(l.xs, l.ys)
        x, _ = to_svg(coord, x, 0)
        if y > 1
            y = 1.0
        elseif y < 0
            y = 0.0
        end
        y = ytop + height * (1 - y)
        x = num2str(x)
        y = num2str(y)
        if @isdefined(last_y) && last_y == y
            is_pending = true
            pending = (x, y)
            continue
        end
        if is_pending
            push!(points, "$(pending[1]),$(pending[2])")
        end
        is_pending = false
        last_y = y
        push!(points, "$x,$y")
    end
    if is_pending
        push!(points, "$(pending[1]),$(pending[2])")
    end
    line["points"] = join(points, " ")
    apply_attrs!(line, l.attrs)
    return
end

function Visual._render_ezxml(io::IO, info::Visual.SVGInfo)
    has_title = !isempty(info.title.title)
    has_plot = !isempty(info.lines)
    doc = get_template(info.trap.name, has_title, has_plot)
    svg = root(doc)
    coord = SVGCoord(svg)
    if has_title
        title = findfirst("//svg:text[@id='Title']", svg, svg_ns)
        @assert title !== nothing
        title.content = info.title.title
        apply_attrs!(title, info.title.attrs)
    end
    for c in info.circles
        circ = addelement!(svg, "circle")
        sx, sy = to_svg(coord, c.x, c.y)
        circ["cx"] = num2str(sx)
        circ["cy"] = num2str(sy)
        circ["r"] = num2str(to_svg(coord, c.r))
        apply_attrs!(circ, c.attrs)
    end
    if has_plot
        plotspace = findfirst("//svg:rect[@id='PlotSpace']", svg, svg_ns)
        @assert plotspace !== nothing
        ytop = parse(Float64, plotspace["y"])
        height = parse(Float64, plotspace["height"])
        unlink!(plotspace)

        for l in info.lines
            draw_line(svg, l, coord, ytop, height)
        end
    end
    if !isempty(info.electrode_colors)
        electrodes = Dict{String,Vector{EzXML.Node}}()
        for ele in findall("//svg:path[@id]", svg, svg_ns)
            push!(get!(()->EzXML.Node[], electrodes, split(ele["id"], "-")[1]), ele)
        end
        for (name, color_str) in info.electrode_colors
            for ele in electrodes[name]
                ele["style"] = replace(ele["style"], r"fill:[^;]*"=>"fill:$color_str")
            end
        end
    end
    println(io, doc)
    return
end

end
