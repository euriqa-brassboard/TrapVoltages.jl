#

module Traps

export TrapDesc

struct TrapDesc
    ele_names::Vector{String}
    ele_indices::Dict{String,Int}
    function TrapDesc(ele_names)
        ele_indices = Dict{String,Int}()
        for (i, name) in enumerate(ele_names)
            ele_indices[name] = i
        end
        return new(ele_names, ele_indices)
    end
end

const _trap_px = TrapDesc(["GND"; "RF";
                           "L" .* string.(0:9);
                           "O" .* string.(0:1);
                           "Q" .* string.(0:65);
                           "S" .* string.(0:11);])
const _trap_hoa = TrapDesc(["GND"; "RF";
                            "G" .* string.(1:8);
                            "L" .* string.(1:16);
                            "Q" .* string.(1:40);
                            "T" .* string.(1:6);
                            "Y" .* string.(1:24);])

TrapDesc(trap::TrapDesc) = trap
function TrapDesc(name::AbstractString)
    if name == "phoenix" || name == "peregrine"
        return _trap_px
    elseif name == "hoa"
        return _trap_hoa
    end
    throw(ArgumentError("Unknown trap name $(name)"))
end

end
