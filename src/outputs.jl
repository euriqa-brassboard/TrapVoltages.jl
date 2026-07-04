#!/usr/bin/julia

"""
Functions for processing the final output files for artiq
"""
module Outputs

import ..Potentials.Potential

export MapFile, CompensationFile, TransferFile, LineMap,
    load_file, write_file, map_line

function with_write(cb, file::IO)
    cb(file)
end

function with_write(cb, file::AbstractString)
    open(cb, file, "w")
end

function str_float(v)
    if isinteger(v)
        return string(Int(v))
    end
    return uppercase(string(v))
end

function eachline_nocomment(file)
    it0 = eachline(file)
    it1 = Iterators.map(line->strip(split(line, "#", limit=2)[1]), it0)
    return Iterators.filter(!isempty, it1)
end

struct MapFile
    # Only the names seem to be useful
    names::Vector{String}
    MapFile(names=String[]) = new(names)
end

function load_file(file, ::Type{MapFile})
    res = MapFile()
    for line in eachline_nocomment(file)
        push!(res.names, split(line, '\t', limit=2)[1])
    end
    return res
end

function write_file(file, data::MapFile)
    with_write(file) do fh
        id = 0
        for name in data.names
            print(fh, "$(name)\t$(id)\t$(id)\r\n")
            id += 1
        end
    end
end

function Base.Dict(mapfile::MapFile)
    res = Dict{String,Int}()
    for i in 1:length(mapfile.names)
        res[mapfile.names[i]] = i
    end
    return res
end

struct LineMap
    electrode_map::Vector{Int}
end

function _assign_index!(electrode_map, potential, name_map, @nospecialize(ele), i)
    if isa(ele, Integer)
        _assign_index!(electrode_map, potential, name_map,
                       potential.electrode_names[ele], i)
    elseif isa(ele, AbstractString)
        electrode_map[name_map[ele]] = i
    elseif isa(ele, AbstractVector)
        for e in ele
            _assign_index!(electrode_map, potential, name_map, e, i)
        end
    else
        throw(TypeError(:LineMap, "Unknown electrode ID type",
                        Union{Integer,AbstractString,AbstractVector}, ele))
    end
end

function LineMap(potential::Potential, mapfile::MapFile, electrodes)
    nelectrodes = length(mapfile.names)
    name_map = Dict(mapfile)
    electrode_map = zeros(Int, nelectrodes)
    for (i, ele) in enumerate(electrodes)
        _assign_index!(electrode_map, potential, name_map, ele, i)
    end
    return LineMap(electrode_map)
end

function map_line(lm::LineMap, term; buff=nothing)
    nelectrodes = length(lm.electrode_map)
    if buff === nothing
        buff = Vector{Float64}(undef, nelectrodes)
    elseif length(buff) != nelectrodes
        throw(ArgumentError("Buffer size mismatch"))
    end
    for i in 1:nelectrodes
        i2 = lm.electrode_map[i]
        if i2 > 0
            buff[i] = term[i2]
        else
            buff[i] = 0
        end
    end
    return buff
end

struct CompensationFile
    map::MapFile
    term_names::Vector{String}
    term_values::Vector{Vector{Float64}}
end

function load_file(file, ::Type{CompensationFile}; mapfile=nothing)
    # Assume the eachline iterator is stateful
    term_names = String[]
    lines = eachline_nocomment(file)
    for line in lines
        fields = split(line, '\t')
        assign_expr = split(fields[1], '=', limit=2)
        if length(assign_expr) < 2
            # We found the last line of the term names
            # and now it's the line with the electrode names
            if mapfile === nothing
                mapfile = MapFile(fields)
            elseif mapfile.names != fields
                throw(ArgumentError("Electrode name mismatch between MapFile and CompensationFile"))
            end
            break
        end
        # Read term names
        push!(term_names, assign_expr[1])
    end
    mapfile = mapfile::MapFile
    # Read the rest of the lines
    term_values = Vector{Float64}[]
    for line in lines
        fields = parse.(Float64, split(line, '\t'))
        if length(fields) != length(mapfile.names)
            throw(ArgumentError("Wrong number of electrode for voltage solution"))
        end
        push!(term_values, fields)
    end
    if length(term_values) != length(term_names)
        throw(ArgumentError("Mismatch between the number of term names and term values"))
    end
    return CompensationFile(mapfile, term_names, term_values)
end

function write_file(file, data::CompensationFile)
    with_write(file) do fh
        names_suffix = '\t' ^ (length(data.map.names) - 1)
        id = 0
        for name in data.term_names
            # \n line end seems to work fine with the compensation file
            print(fh, "$(name)=$(id)$(names_suffix)\n")
            id += 1
        end
        println(fh, join(data.map.names, '\t'))
        for values in data.term_values
            println(fh, join((str_float(val) for val in values), '\t'))
        end
    end
end

struct TransferFile
    map::MapFile
    line_values::Vector{Vector{Float64}}
end

function load_file(file, ::Type{TransferFile}; mapfile=nothing)
    # Assume the eachline iterator is stateful
    lines = eachline_nocomment(file)
    names = split(first(lines), '\t')
    if mapfile === nothing
        mapfile = MapFile(names)
    elseif mapfile.names != names
        throw(ArgumentError("Electrode name mismatch between MapFile and TransferFile"))
    end
    mapfile = mapfile::MapFile
    # Read the rest of the lines
    line_values = Vector{Float64}[]
    for line in lines
        fields = parse.(Float64, split(line, '\t'))
        if length(fields) != length(mapfile.names)
            throw(ArgumentError("Wrong number of electrode for voltage solution"))
        end
        push!(line_values, fields)
    end
    return TransferFile(mapfile, line_values)
end

function write_file(file, data::TransferFile)
    with_write(file) do fh
        print(fh, join(data.map.names, '\t') * "\r\n")
        for values in data.line_values
            # The original file lacks end-of-file new line
            # Hopefully that doesn't matter as much.
            print(fh, join((str_float(val) for val in values), '\t') * "\r\n")
        end
    end
end

end
