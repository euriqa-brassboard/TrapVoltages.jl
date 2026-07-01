module TrapVoltages

function gradient end
function get_single end

include("poly_fit.jl")

include("traps.jl")
using .Traps
export TrapDesc
public find_electrodes

include("units.jl")

include("potentials.jl")
include("outputs.jl")

include("optimizers.jl")

include("solutions.jl")

end # module TrapVoltages
