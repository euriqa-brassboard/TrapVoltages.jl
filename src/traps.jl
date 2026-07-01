#

module Traps

export TrapDesc

struct ElectrodePosition
    name::String
    left::Float64
    right::Float64
    up::Bool
end

struct RegionElectrodePositions
    inner::Vector{ElectrodePosition}
    outer::Vector{ElectrodePosition}
end

struct TrapDesc
    name::String
    ele_names::Vector{String}
    ele_indices::Dict{String,Int}
    ele_region_pos::Vector{RegionElectrodePositions}
    function TrapDesc(name, ele_names)
        ele_indices = Dict{String,Int}()
        for (i, ele_name) in enumerate(ele_names)
            ele_indices[ele_name] = i
        end
        return new(name, ele_names, ele_indices, RegionElectrodePositions[])
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

function __populate_phoenix_positions(trap)
    # X position of electrodes for phoenix and peregrine traps.
    # We need to know the axial (X) positions of the electrodes so that we can figure out
    # which electrodes to use for generating potentials at a given location.

    # From looking at the Sandia 3D model, each inner electrode is 70um wide in X direction
    # (67um + 3um gap) and each outer quantum is 2x this (140um total).
    # All the odd electrode are always located at the same X position as
    # the electrode numbered one less than it.

    # In unit of 70um and showing only the even electrodes, the order/positions
    # of the electrodes are,

    # Outer: |               46.5(O0)            |22(Q44-64)|          14.5(O0)         |
    # Inner: |2(gap)|10(GND)|5(L0-8)|30(S10-0 x 5)|22(Q0-42) |8(S0-10,0-3)|6(GND)|2(gap)|

    # where the number outside the parenthesis is the width in unit of 70um
    # and the parenthesis marks the corresponding (even) electrode.
    # The origin is located in the middle of the quantum region (11 from left and right).
    # This distribution is cross-checked with the potential data
    # by setting two pairs of electrode to opposite values and measuring the position
    # of the zero crossing.

    begin_gnd = 12
    nL = 5
    nS = 6
    S_rep1 = 5
    nQ = 22
    S_rep2 = 1
    end_gnd = 8

    unit_um = 70

    @assert nQ % 2 == 0
    nQ_outer = nQ ÷ 2
    left_edge = -(begin_gnd + nL + nS * S_rep1 + nQ ÷ 2)

    pos_inner = left_edge + begin_gnd
    pos_outer = left_edge

    positions = RegionElectrodePositions(ElectrodePosition[], ElectrodePosition[])
    push!(trap.ele_region_pos, positions)

    # Loading
    for i in 1:nL
        push!(positions.inner,
              ElectrodePosition("L$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(positions.inner,
              ElectrodePosition("L$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Transition 1
    for j in 1:S_rep1
        for i in nS:-1:1
            push!(positions.inner,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, true))
            push!(positions.inner,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, false))
            pos_inner += 1
        end
    end

    # Outer 1
    push!(positions.outer,
          ElectrodePosition("O0", pos_outer * unit_um,
                            (pos_inner - 0.5) * unit_um, true))
    push!(positions.outer,
          ElectrodePosition("O1", pos_outer * unit_um,
                            (pos_inner - 0.5) * unit_um, false))
    pos_outer = pos_inner - 0.5

    # Quantum inner
    for i in 1:nQ
        push!(positions.inner,
              ElectrodePosition("Q$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(positions.inner,
              ElectrodePosition("Q$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Quantum outer
    for i in 1:nQ_outer
        i += nQ
        push!(positions.outer,
              ElectrodePosition("Q$(i * 2 - 2)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um, true))
        push!(positions.outer,
              ElectrodePosition("Q$(i * 2 - 1)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um, false))
        pos_outer += 2
    end
    @assert pos_inner - 0.5 == pos_outer

    # Transition 2
    for j in 1:S_rep2
        for i in 1:nS
            push!(positions.inner,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, true))
            push!(positions.inner,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, false))
            pos_inner += 1
        end
    end

    # S0-S3 appeared again at the end (shouldn't really matter......)
    for i in 1:4
        push!(positions.inner,
              ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(positions.inner,
              ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Outer 2
    push!(positions.outer,
          ElectrodePosition("O0", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um, true))
    push!(positions.outer,
          ElectrodePosition("O1", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um, false))
    return
end
__populate_phoenix_positions(_trap_phoenix)
__populate_phoenix_positions(_trap_peregrine)

end
