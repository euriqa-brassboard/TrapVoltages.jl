module TrapVoltages

function gradient end
function get_single end

include("poly_fit.jl")

include("traps.jl")
using .Traps
export TrapDesc

include("potentials.jl")
include("outputs.jl")

end # module TrapVoltages
