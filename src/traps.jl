#

module Traps

export TrapDesc

struct TrapDesc
    name::String
    ele_names::Vector{String}
    ele_indices::Dict{String,Int}
    function TrapDesc(name, ele_names)
        ele_indices = Dict{String,Int}()
        for (i, ele_name) in enumerate(ele_names)
            ele_indices[ele_name] = i
        end
        return new(name, ele_names, ele_indices)
    end
end

const _trap_phoenix = TrapDesc("phoenix", ["GND"; "RF";
                                "L" .* string.(0:9);
                                "O" .* string.(0:1);
                                "Q" .* string.(0:65);
                                "S" .* string.(0:11);])
const _trap_peregrine = TrapDesc("peregrine", ["GND"; "RF";
                                  "L" .* string.(0:9);
                                  "O" .* string.(0:1);
                                  "Q" .* string.(0:65);
                                  "S" .* string.(0:11);])
const _trap_hoa = TrapDesc("hoa", ["GND"; "RF";
                            "G" .* string.(1:8);
                            "L" .* string.(1:16);
                            "Q" .* string.(1:40);
                            "T" .* string.(1:6);
                            "Y" .* string.(1:24);])

TrapDesc(trap::TrapDesc) = trap
function TrapDesc(name::AbstractString)
    if name == "phoenix"
        return _trap_phoenix
    elseif name == "peregrine"
        return _trap_peregrine
    elseif name == "hoa"
        return _trap_hoa
    end
    throw(ArgumentError("Unknown trap name $(name)"))
end

end
