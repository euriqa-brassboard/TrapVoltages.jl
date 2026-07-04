module TrapVoltages

include("utils.jl")

function gradient end
function get_single end

include("poly_fit.jl")

include("traps.jl")
using .Traps: TrapDesc, find_electrodes
export TrapDesc
public find_electrodes

include("units.jl")
using .Units: TrapUnits, MHz_unit
public TrapUnits, MHz_unit

include("potentials.jl")
using .Potentials: import_pillbox_v0, import_pillbox_v1, import_pillbox_64,
    Potential, get_potential, Fitting, get_electrodes, load_short_map,
    x_index_to_axis, y_index_to_axis, z_index_to_axis,
    x_axis_to_index, y_axis_to_index, z_axis_to_index
export import_pillbox_v0, import_pillbox_v1, import_pillbox_64,
    get_potential, load_short_map
public Potential, Fitting, get_electrodes,
    x_index_to_axis, y_index_to_axis, z_index_to_axis,
    x_axis_to_index, y_axis_to_index, z_axis_to_index

include("outputs.jl")

include("optimizers.jl")
using .Optimizers: optimize_minmax, optimize_minmax_span
export optimize_minmax, optimize_minmax_span

include("solutions.jl")
using .Solutions: find_flat_point, find_all_flat_points, CenterTracker, TermMask,
    compensate_terms, solve_compensate, solve_target
public find_flat_point, find_all_flat_point, CenterTracker, TermMask,
    compensate_terms, solve_compensate, solve_target

include("visual.jl")

end # module TrapVoltages
